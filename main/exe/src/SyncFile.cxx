//
// Created by paschen on 3/7/19.
//

#include <iostream>
#include <string>
#include <queue>
#include "eudaq/OptionParser.hh"
#include "eudaq/FileReader.hh"
#include "eudaq/FileWriter.hh"
#include "eudaq/PluginManager.hh"

#include "TFile.h"
#include "TTree.h"


class SyncSetup {
public:
    // event queues to work as independent FIFOs (First In First Out)
    std::queue<eudaq::DetectorEvent> events;  // holds global events
    typedef std::queue<std::shared_ptr<eudaq::Event>> eventqueue_t;  // queue holding shared pointers to eudaq events
    std::vector<eventqueue_t> subevents;  // vector of subevent queues

    std::vector<unsigned> triggerIDs;
    std::vector<unsigned> triggerID_offsets;
    std::vector<unsigned> typeIDs;
    std::vector<std::string> subtype;

    std::vector<bool> high_trg_nr_problem;
    unsigned max_triggerID = 0;
    size_t n_sub;

    unsigned tlu_id = eudaq::Event::str2id("_TLU");  // ID of TLU subevents

    // constructor needs BORE (Begin Of Run Event) to set up the queues
    explicit SyncSetup(const eudaq::DetectorEvent *bore) {
        if (!bore->IsBORE())
            throw;
        // determine number of sub detectors/events
        n_sub = bore->NumEvents();
        // initialize all counters and subevents
        triggerIDs.resize(n_sub, 0);
        triggerID_offsets.resize(n_sub, 0);
        high_trg_nr_problem.resize(n_sub, false);
        subevents.resize(n_sub);
        for (size_t i = 0; i < n_sub; i++) {
            typeIDs.push_back(bore->GetEvent(i)->get_id());
            subtype.push_back(bore->GetEvent(i)->GetSubType());
        }
    }

    // checks for TLU trigger number overflow/wraparound (15 bit -> occurs every 32768 triggers)
    void check_tlu_overflow(size_t i_subevent, unsigned prev_trigger) {
        size_t i = i_subevent;
        // TODO: define detection limit, currently trying to detect jumps of up to ~500-1000
        if (triggerIDs[i] < std::max((long) triggerID_offsets[i] + 500, (long) prev_trigger - 32000) &&
                prev_trigger > triggerIDs[i] + 32000) {
            std::cout << "overflow detected " << subtype[i] << " new: " << triggerIDs[i] << " old: "
                      << prev_trigger << std::endl;
            triggerID_offsets[i] += 32768;
            triggerIDs[i] = triggerID_offsets[i] + eudaq::PluginManager::GetTriggerID(*subevents[i].front());
        }
    }

    // updates trigger ID
    // checks for FEI4 - high trigger number problem (>32768)
    // calls trigger number overflow/wraparound check
    void update_trigger_id(size_t i) {
        unsigned triggerID_old = triggerIDs[i];
        // check for FEI4 problem, do not increase trigger to not loose synchronization.
        // most triggers seem to be sent later (checked for DEPFET 2019 run 002141)
        unsigned cur_trigID = eudaq::PluginManager::GetTriggerID(*subevents[i].front());
        if (typeIDs[i] != tlu_id && cur_trigID > 32768) {
            std::cout << eudaq::Event::id2str(typeIDs[i]) << " " << subtype[i] << " " << cur_trigID << std::endl;
            // triggerIDs[i]++;
            high_trg_nr_problem[i] = true;

        // catch NI error and increase normally
//        } else if (subtype[i] == "NI" && (cur_trigID == 4800 || std::abs((long) (cur_trigID + triggerID_offsets[i]) - (long)triggerID_old) > 5)) {
//            std::cout << "NI problem " << cur_trigID + triggerID_offsets[i] << " " << triggerID_old <<  std::endl;
//            triggerIDs[i]++;
//            high_trg_nr_problem[i] = true;
        } else {
            // reset trg_nr_problem flag
            if (high_trg_nr_problem[i])
                high_trg_nr_problem[i] = false;
            triggerIDs[i] = triggerID_offsets[i] + cur_trigID;
            // check for overflow
            check_tlu_overflow(i, triggerID_old);
        }
    }

    void add_event(const eudaq::DetectorEvent *read_event) {

        // get run number, event number and timestamp for copying
        unsigned runnumber = read_event->GetRunNumber();
        unsigned eventnumber = read_event->GetEventNumber();
        uint64_t timestamp = read_event->GetTimestamp();

        // create global event in queue
        events.emplace(eudaq::DetectorEvent(runnumber, eventnumber, timestamp));
        // copy flags
        events.back().SetFlags(read_event->GetFlags());
        // copy tags
        std::vector<std::string> taglist = read_event->GetTagList("");
        for (const auto & i : taglist) {
            events.back().SetTag(i, read_event->GetTag(i));
        }
        // read in sub events
        for (size_t i = 0; i < n_sub; i++) {
            // create shared pointer references to subevents in queues
            subevents[i].emplace(std::shared_ptr<eudaq::Event>(read_event->GetEventPtr(i)));
            // if queue was empty, determine new trigger ID
            if (subevents[i].size() == 1) {
                update_trigger_id(i);
            }
        }
    }

    // checks for FEI4 high trigger number problem
    bool check_high_trg() {
        bool any_high_trg_nr_problem = false;
        for (size_t i = 0; i < n_sub; i++)
            any_high_trg_nr_problem |= high_trg_nr_problem[i];
        if (any_high_trg_nr_problem) {
            std::cout << events.front() << std::endl;
            events.pop();
            for (size_t i = 0; i < n_sub; i++) {
                std::cout << *subevents[i].front() << std::endl;
                subevents[i].pop();
            }
            std::cout << "High trigger number problem, skip event" << std::endl << std::endl;
            return true;
        }
        return false;
    }

    void update_max_trigger() {
        // TODO: define max trigger jump
        for (size_t i = 0; i < n_sub; i++) {
            if ((typeIDs[i] != tlu_id) && (triggerIDs[i] > max_triggerID) && (triggerIDs[i] < max_triggerID + 10000) && !high_trg_nr_problem[i]) {
                max_triggerID = triggerIDs[i];
                // mismatched_det = i;
            }
        }
    }

    // checks event for mismatch
    bool check_event() {
        // check trigger IDs for mismatch
        for (size_t i = 0; i < n_sub; i++) {
            if (typeIDs[i] != tlu_id && (triggerIDs[i] != events.front().GetEventNumber() + 1 || high_trg_nr_problem[i])) {  // && !subevents[i].front()->IsEORE()) {
                std::cout << i << " " << eudaq::Event::id2str(typeIDs[i]) << ":" << subtype[i] << " "
                          << triggerIDs[i] << std::endl;
                return true;
            }
        }
        return false;
    }

    // discard global events and subevents of lower trigger IDs
    void discard_mismatch() {
        std::cout << "max_trigger_ID " << max_triggerID << std::endl;
        std::cout << events.front() << std::endl;
        for (size_t i = 0; i < n_sub; i++)
            std::cout << *subevents[i].front() << std::endl;
        // check global event
        if ((events.front().GetEventNumber() + 1) != max_triggerID) events.pop();
        // check sub events
        for (size_t i = 0; i < n_sub; i++) {
            // check tlu event
            if (typeIDs[i] == tlu_id) {
                if ((subevents[i].front()->GetEventNumber() + 1) != max_triggerID) {
                    std::cout << "delete TLU subevent " << subevents[i].front()->GetEventNumber() << std::endl;
                    subevents[i].pop();
                }
            }
            // check other sub events
            else if (triggerIDs[i] != max_triggerID || high_trg_nr_problem[i]) {
                std::cout << "delete " << i << " " << eudaq::Event::id2str(typeIDs[i]) << ":" << subtype[i]
                          << " " << triggerIDs[i] << std::endl;
                subevents[i].pop();
                if (!subevents[i].empty()) {
                    update_trigger_id(i);
                }
            }
        }
    }

    // merges current subevents into global event and deletes subevents from their respective queues
    eudaq::DetectorEvent &merge_event() {
        for (size_t i = 0; i < n_sub; i++) {
            events.front().AddEvent(subevents[i].front());
            // subevent_trigger[i] = eudaq::PluginManager::GetTriggerID(*subevents[i].front());
            // delete first event of subevent queue
            subevents[i].pop();
            if (!subevents[i].empty()) {
                update_trigger_id(i);
            }
        }
        return events.front();
    }

    void pop_event() {
        // write out global event in deque and delete
        events.pop();
    }
} ;


int main(int argc, char **argv) {
    eudaq::OptionParser op("EUDAQ RAW data syncer", "1.0",
                           "A command-line tool for syncing sub events based on TLU trigger numbers", -1);
    eudaq::Option<std::string> infile(op, "i", "infile", "", "string", "Input filename");
    eudaq::Option<std::string> outpath(op, "o", "outpath", "", "string", "Output path");
    try {

        op.Parse(argv);
        // check options
        if (infile.Value().empty()) {
            throw (eudaq::OptionException("Specify input file"));
        }
        if (outpath.Value().empty()) {
            throw (eudaq::OptionException("Specify output path"));
        }
        // remove trailing slashes from path name
        std::string output_file = outpath.Value();
        while (output_file.back() == '/' || output_file.back() == '\\') {
            output_file.erase(output_file.size() - 1);
        }
        // get file name of input file
        size_t i = infile.Value().rfind('/', infile.Value().length());
        if (i == std::string::npos) {
            i = infile.Value().rfind('\\', infile.Value().length());
        }
        // FIXME: only works for linux so far
        if (i != std::string::npos)
            output_file = output_file + '/' + infile.Value().substr(i + 1, infile.Value().length() - i);
        else
            output_file = output_file + '/' + infile.Value();

        // reader instance opening file
        eudaq::FileReader reader(infile.Value());
        std::cout << "Reading file: " << reader.Filename() << std::endl;  // display file name
        // pointer to constant detector event
        eudaq::DetectorEvent const *read_event;

        // writer instance, opening file
        eudaq::FileWriter *writer = eudaq::FileWriterFactory::Create("native");
        writer->SetFilePattern(output_file);
        writer->StartRun(reader.RunNumber());
        std::cout << "Writing to file: " << output_file << std::endl;

        // open root file and define tree
        TFile rootfile((output_file.substr(0, output_file.length()-3)+"root").c_str(), "RECREATE");
        std::cout << "Writing root file: " << rootfile.GetName() << std::endl;
        // Initialize tree
        TTree event_tree("events", "Events");
        unsigned eventnumber;
        event_tree.Branch("Event", &eventnumber, "Event/i");

        // get first event (BORE - Beginning Of Run Event) and copy to output file
        read_event = &reader.GetDetectorEvent();
        writer->WriteEvent(*read_event);
        // initialize PluginManager with BORE
        eudaq::PluginManager::Initialize(*read_event);
        // initialize SyncSetup
        SyncSetup setup(read_event);

        // variables used in the loop
        unsigned counter = 0;
        unsigned end = 1000000;

        unsigned subevent_trigger[setup.n_sub];
        // for (size_t i = 0; i < n_sub; i++) event_tree.Branch("subevent_trigger", &subevent_trigger[i], "subevent_trigger/i");
        // event_tree.Branch("subevent_trigger", subevent_trigger, ("subevent_trigger["+std::to_string(setup.n_sub)+"]/i").c_str());


        // mismatch variables
        bool got_mismatch = false;
        size_t mismatched_det = 0;

        while (reader.NextEvent()) { // && counter < 37000) {
            counter++;
            // get current event as detector event
            read_event = &reader.GetDetectorEvent();
            // TODO handle EORE, also in subevents...
//            if(read_event->IsEORE()) {
//                writer->WriteEvent(*read_event);
//                continue;
//            }

            // fill global + subevent queues and determine total trigger IDs
            setup.add_event(read_event);

            // if FEI4 trigger number problem, immediately discard entire event
//            if(setup.check_high_trg())
//                 continue;

            // check event for mismatch
            got_mismatch = setup.check_event();
            setup.update_max_trigger();

            // in case of mismatch discard global events and subevents of lower trigger IDs
            if (got_mismatch) {
                std::cout << "reset mismatch" << std::endl;
                setup.discard_mismatch();
            }
            // no mismatch: write out event
            else {
                writer->WriteEvent(setup.merge_event());
                setup.pop_event();
            }
        }

        // TODO: Write out possibly remaining events

        std::cout << "Remaining global events: " << setup.events.size() << std::endl;
        for (size_t i = 0; i < setup.n_sub; i++) {
            std::cout << "Remaining subevents: " << std::endl;
            std::cout << eudaq::Event::id2str(setup.typeIDs[i]) << " " << setup.subtype[i] << " " << setup.subevents[i].size()
                      << std::endl;
        }
    event_tree.Write();
    rootfile.Close();

    } catch (...) {
        return op.HandleMainException();
    }

    return 0;
}
