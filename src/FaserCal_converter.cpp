#include <fstream>
#include <vector>
#include <iostream>
#include <cstring>
#include "EventFormats/OCBDecoder.hpp"
#include <TFile.h>
#include <TTree.h>

static uint32_t bytes_to_uint32(const unsigned char buf[4]) {
    return (static_cast<uint32_t>(buf[0]) )
        | (static_cast<uint32_t>(buf[1]) << 8)
        | (static_cast<uint32_t>(buf[2]) << 16)
        | (static_cast<uint32_t>(buf[3]) << 24);
}

struct FEBhit {
    int m_channel_id;
    int m_amplitude_lg;
    int m_amplitude_hg;
    std::vector<int> m_hit_time_rise;
    std::vector<int> m_hit_time_fall;
    std::vector<int> m_gts_time_rise;
    std::vector<int> m_gts_time_fall;

    FEBhit(int id) : m_channel_id(id) {};
};

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <binary-file> <output-file> [n_max] [debug] [new_fw]\n";
        return 1;
    }

    const char* path = argv[1];
    const char* path_output = argv[2];

    bool new_fw = true;
    bool debug = false;
    int n_max = -1;

    if (argc > 3) {
        std::string argv_2 = argv[3];
        n_max = std::stoi(argv_2);
    }
    if (argc > 4) {
        std::string argv_3 = argv[4];
        debug = (argv_3 == "true" || argv_3 == "1");
    }
    if (argc > 5) {
        std::string argv_4 = argv[5];
        new_fw = (argv_4 == "true" || argv_4 == "1");
    }

    std::ifstream in(path, std::ios::binary);
    if (!in) {
        std::cerr << "Failed to open file: " << path << "\n";
        return 2;
    }

    int index = 0;
    int start_index = -1;
    bool is_open_packet = false;
    std::vector<uint32_t> word_list;

    int glob_evt_counter = 0;
    int glob_hit_counter = 0;
    int glob_time_hit_counter = 0;
    int glob_empty_packet_counter = 0;

    TFile f(path_output, "RECREATE");
    TTree tree("Events", "FASERCal events");

    // Define branches for TTree

    // single values per event
    int bcid;
    int event_nbr;
    int event_nbr_extended;

    // arrays of single values, 1 per hit (amplitude-like info)
    std::vector<int> channel_id;
    std::vector<int> feb_id;
    std::vector<int> amplitude_lg;
    std::vector<int> amplitude_hg;

    // arrays of arrays, 1 array of values per hit (time-like info: allow more than 1 time hit per board and channel)
    std::vector<std::vector<int>> hit_time_rise;
    std::vector<std::vector<int>> hit_time_fall;
    std::vector<std::vector<int>> gts_time_rise;
    std::vector<std::vector<int>> gts_time_fall;
    tree.Branch("bcid", &bcid);
    tree.Branch("event_nbr", &event_nbr);
    tree.Branch("event_nbr_extended", &event_nbr_extended);
    tree.Branch("feb_id", &feb_id);
    tree.Branch("channel_id", &channel_id);
    tree.Branch("amplitude_lg", &amplitude_lg);
    tree.Branch("amplitude_hg", &amplitude_hg);
    tree.Branch("hit_time_rise", &hit_time_rise);
    tree.Branch("hit_time_fall", &hit_time_fall);
    tree.Branch("gts_time_rise", &gts_time_rise);
    tree.Branch("gts_time_fall", &gts_time_fall);

    unsigned char buf[4];
    while (in.read(reinterpret_cast<char*>(buf), 4)) {
        uint32_t w = bytes_to_uint32(buf);
        word_list.push_back(w);
        std::unique_ptr<Word> parsed = parse_word(w);
        std::cout << *parsed;
        if (parsed->word_id == WordID::OCB_PACKET_HEADER && !is_open_packet) {
            start_index = index;
            is_open_packet = true; // This allows first OCB header to also be added to the OCB packet in the new FW
        }
        if (parsed->word_id == WordID::OCB_PACKET_TRAILER) {
            if (start_index < 0) {
                throw std::runtime_error("OCB Packet Trailer received without corresponding Header");
            }

            std::vector<uint32_t> ocb_packet_word_list;
            for (int k = start_index; k <= index; ++k) {
                ocb_packet_word_list.push_back(word_list[static_cast<size_t>(k)]);
            }

            // construct OCBDataPacket from pointer and size in bytes
            OCBDataPacket ev = OCBDataPacket(ocb_packet_word_list.data(), ocb_packet_word_list.size() * 4, debug, new_fw);
            std::cout << ev;

            // Clear all vectors used to fill in the TTree
            channel_id.clear();
            feb_id.clear();
            amplitude_lg.clear();
            amplitude_hg.clear();
            hit_time_rise.clear();
            hit_time_fall.clear();
            gts_time_rise.clear();
            gts_time_fall.clear();

            // Retrieve event information for ntuples
            event_nbr = ev.get_event_id();
            bcid = ev.get_bcid();
            event_nbr_extended = ev.get_event_id_extended();

            int n_feb_with_data = 0;

            for (size_t board_id = 0; board_id < OCBConfig::NUM_FEBS_PER_OCB; board_id++){
                if (ev.isCorruptedFEB(board_id)) {
                    std::cout << "FEB " << board_id << " data packet is corrupted (missing header or trailer)." << std::endl;
                }
                if (ev.hasData(board_id)) {
                    std::cout << "FEB " << board_id << " sent data." << std::endl;

                    auto feb_packet = ev[board_id];
                    if (feb_packet.get_hit_amplitudes().size() > 0) n_feb_with_data++;
                    if (feb_packet.isCorrupted()) { // this should never happen if the FEB is different from nullptr because the check is also made at the OCB level
                        std::cout << "FEB " << board_id << " data packet is corrupted (missing header or trailer)." << std::endl; 
                    }
                    if (feb_packet.hasMissingGTS()) {
                        std::cout << "FEB " << board_id << " data packet has missing GTS header or trailer." << std::endl;
                    }

                    // Create one Hit struct for each fired channel
                    std::map<int, FEBhit> feb_hit_map;
                    for (const auto& hit_amplitude : feb_packet.get_hit_amplitudes()) {
                        int hit_channel_id = hit_amplitude.get_channel_id();
                        if (feb_hit_map.find(hit_channel_id) == feb_hit_map.end()) {
                            feb_hit_map.emplace(hit_channel_id, FEBhit(hit_channel_id));
                            auto it = feb_hit_map.find(hit_channel_id);
                            it->second.m_amplitude_lg = hit_amplitude.get_amplitude_lg();
                            it->second.m_amplitude_hg = hit_amplitude.get_amplitude_hg();
                            glob_hit_counter++;
                        }
                        else {
                            throw std::runtime_error("Duplicate hit for channel " + std::to_string(hit_channel_id) + " in FEB " + std::to_string(board_id));
                        }
                    }

                    for (const auto& hit_time : feb_packet.get_hit_times()) {
                        int hit_channel_id = hit_time.get_channel_id();
                        if (feb_hit_map.find(hit_channel_id) != feb_hit_map.end()) { 
                            // NOTE: if a time hit is detected without corresponding amplitude value for the same channel, such hit is not saved
                            auto it = feb_hit_map.find(hit_channel_id);
                            it->second.m_hit_time_rise.push_back(hit_time.get_hit_time_rise());
                            it->second.m_hit_time_fall.push_back(hit_time.get_hit_time_fall());
                            it->second.m_gts_time_rise.push_back(hit_time.get_gts_tag_rise());
                            it->second.m_gts_time_fall.push_back(hit_time.get_gts_tag_fall());
                            glob_time_hit_counter++;
                        }
                    }

                    for (const auto& hit_pair : feb_hit_map) {
                        const FEBhit& hit = hit_pair.second;
                        // Fill TTree branches with hit data
                        channel_id.push_back(hit.m_channel_id);
                        feb_id.push_back(board_id);
                        amplitude_lg.push_back(hit.m_amplitude_lg);
                        amplitude_hg.push_back(hit.m_amplitude_hg);
                        hit_time_rise.push_back(hit.m_hit_time_rise);
                        hit_time_fall.push_back(hit.m_hit_time_fall);
                        gts_time_rise.push_back(hit.m_gts_time_rise);
                        gts_time_fall.push_back(hit.m_gts_time_fall);
                    }
                }
            }


            tree.Fill();

            if (n_feb_with_data == 0) {
                glob_empty_packet_counter++;
            }

            start_index = -1;
            is_open_packet = false;
            glob_evt_counter++;
        }

        if (glob_evt_counter >= n_max && n_max > 0) {
            break;
        }
        ++index;
    }
    std::cout << "Number of OCB packets: " << glob_evt_counter << std::endl;
    std::cout << "Number of empty OCB packets: " << glob_empty_packet_counter << std::endl;
    std::cout << "Number of hits: " << glob_hit_counter << std::endl;
    std::cout << "Number of time hits: " << glob_time_hit_counter << std::endl;

    tree.Write();

    return 0;
}
