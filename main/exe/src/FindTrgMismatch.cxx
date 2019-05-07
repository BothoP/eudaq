//
// Created by paschen on 11/5/18.
//

#include <iostream>
#include <fstream>
#include <stdio.h>
#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"


int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "usage: " << argv[0] << " file" << std::endl;
        return 1;
    }
    // check if file has size != 0
    std::ifstream in(argv[1], std::ifstream::ate | std::ifstream::binary);
    // std::cout << in.tellg() << std::endl;
    if (in.tellg() == 0) {
        return 1;
    }
    in.close();

    eudaq::FileReader reader(argv[1]);

    // pointer to constant detector event
    eudaq::DetectorEvent const *det_event;
    eudaq::Event const *event;
    eudaq::StandardPlane std_plane;

    // display file name
    std::cout << "Opened file: " << reader.Filename() << std::endl;
    // get first event (BORE - Beginning Of Run Event)
    det_event = &reader.GetDetectorEvent();
    // initialize PluginManager with BORE
    eudaq::PluginManager::Initialize(*det_event);
    // get runnumber from filename
    int substr_index = reader.Filename().rfind("run");
    std::string runnumber_str_from_filename = reader.Filename().substr(substr_index+3, 6);
    std::cout << "Run Number: " << runnumber_str_from_filename << std::endl;

    int runnumber = det_event->GetRunNumber();
    int counter = 0;
    int tlu_event = 0;
    unsigned tlu_id = eudaq::Event::str2id("_TLU");  // ID of TLU subevents
    bool got_mismatch[6] = {false, false, false, false, false, false};
    int trg_mismatch[6] = {-1, -1, -1, -1, -1, -1};
    std::string type[6];

    // run number mismatch?
    if (runnumber != std::stoi(runnumber_str_from_filename)) {
        type[5] = "runnumber";
        got_mismatch[5] = true;
        // trg_mismatch[5] = -1;
    }

    // loop over the file now
    while (reader.NextEvent()) {
        // get current event as detector event
        det_event = &reader.GetDetectorEvent();
        if (det_event->IsEORE()) continue;
        counter++;
        // get current event without casts
        // event = &reader.GetEvent();
        /* Event type << det_event->GetType()      -> yields "DetectorEvent"
           det_event->id2str(det_event->get_id())  -> yields "_DET"
           det_event->GetRunNumber()               -> run number
           det_event->GetEventNumber()             -> event number
           det_event->NumEvents()                  -> number of sub events in detector event
        */

        // get number of event (should be same as TLU trigger number)
        tlu_event = det_event->GetEventNumber();

        // loop over sub events
        for (size_t i = 0; i < det_event->NumEvents(); i++) {
            // get sub event
            eudaq::Event const *sub_event = det_event->GetEvent(i);
            // get trigger ID of sub event
            unsigned triggerID = eudaq::PluginManager::GetTriggerID(*sub_event);
            // get type id of sub event
            unsigned id = sub_event->get_id();

            // check if there is a mismatch
            // sub detectors only take 15 bit trigger numbers and are always ahead of "event number"/tlu trigger by one
            // if first mismatch is detected it's written into the variables
            if ((id != tlu_id) && (triggerID != (tlu_event + 1) % 32768) && !got_mismatch[i]) {
                got_mismatch[i] = true;
                trg_mismatch[i] = tlu_event;
                type[i] = sub_event->GetSubType();
                // std::cout << *det_event << std::endl;
            }
        }
        // check if runnumber is still correct
        runnumber = det_event->GetRunNumber();
        if (runnumber != std::stoi(runnumber_str_from_filename) && !got_mismatch[5]) {
            type[5] = "runnumber";
            got_mismatch[5] = true;
            trg_mismatch[5] = tlu_event;
        }
        if (!(tlu_event % 50000)) {
            std::cout << tlu_event << std::endl;
        }

    }

    // print out mismatches
    for (int i : trg_mismatch) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
    // print event counter
    std::cout << counter << std::endl;

    // open csv file
    std::ofstream csvfile;
    char s[10];
    sprintf(s, "%06d", runnumber);
    std::string srunnumber = s;
    //FIXME: I changed this line to the runnumber from the filesname
    csvfile.open("/run/media/paschen/Seagate Backup Plus Drive/TB2019/eudaq/stats2/mismatch_" + runnumber_str_from_filename + ".csv");

    // get list of indices for detectors with mismatches
    std::vector<int> index;
    for (int i=0; i<6 ; i++) {
        if (got_mismatch[i]) {
            index.push_back(i);
        }
    }

    // write first mismatched detector events to file
    for (auto i: index) {
        csvfile << type[i] << ",";
    }
    csvfile << "number of events\n";
    for (auto i: index) {
        csvfile << trg_mismatch[i] << ",";
    }

    // write number of events in file
    csvfile << counter;
    csvfile.close();
    return 0;
}
