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
    // new events are added to the back
    // current events are read at the front and removed (queue.pop()) once they are matched or discarded
    std::queue<eudaq::DetectorEvent> events;  // holds global events for storing matching sub events to file
    typedef std::queue<std::shared_ptr<eudaq::Event>> eventqueue_t;  // queue holding shared pointers to eudaq events
    std::vector<eventqueue_t> subevents;  // vector of subevent queues

    std::vector<unsigned> trgnr_actual;
    std::vector<unsigned> trgnr_tlu_offset;
    std::vector<unsigned> typeID;
    std::vector<std::string> subtype;

    std::vector<bool> trgnr_wraparound;
    std::vector<bool> trgnr_out_of_order;
    unsigned trgnr_actual_global = 0;
    size_t n_sub;

    unsigned tlu_id = eudaq::Event::str2id("_TLU");  // ID of TLU subevents

    // constructor needs BORE (Begin Of Run Event) to set up the queues
    explicit SyncSetup(const eudaq::DetectorEvent *bore) {
        if (!bore->IsBORE())
            throw;
        // determine number of sub detectors/events
        n_sub = bore->NumEvents();
        // initialize all counters and subevents
        trgnr_actual.resize(n_sub, 0);
        trgnr_tlu_offset.resize(n_sub, 0);
        trgnr_out_of_order.resize(n_sub, false);
        trgnr_wraparound.resize(n_sub, false);
        subevents.resize(n_sub);
        for (size_t i = 0; i < n_sub; i++) {
            typeID.push_back(bore->GetEvent(i)->get_id());
            subtype.push_back(bore->GetEvent(i)->GetSubType());
        }
    }

    bool all_queues_filled() {
        bool any_subqueue_empty = false;
        for (auto &subevent : subevents)
            any_subqueue_empty |= subevent.empty();
        return !(events.empty() || any_subqueue_empty);
    }

    // checks for TLU trigger number overflow/wraparound (15 bit -> occurs every 32768 triggers)
    void check_tlu_overflow(size_t i_subevent, unsigned prev_trigger) {
        size_t i = i_subevent;
        // TODO: define detection limit, currently detecting all backward jumps greater allowed_trg_jump
        // FEI4 trigger jumps have to be detected and corrected before
        // This limit needs to be at least the same size as the NI mismatch detection limit
        static const unsigned allowed_trg_jump = 20;
        if (trgnr_actual[i] + allowed_trg_jump < prev_trigger) {
            std::cout << "overflow detected " << subtype[i] << " new: " << trgnr_actual[i] << " old: "
                      << prev_trigger << std::endl;
            trgnr_tlu_offset[i] += 32768;
            trgnr_actual[i] += 32768;
        }
    }

    // updates trigger ID
    // checks for FEI4 - high trigger number problem (>32768)
    // calls trigger number overflow/wraparound check
    void update_trigger_id(size_t i) {
        // detection limit for NI trigger jumps
        static const unsigned allowed_trg_jump_NI = 5;
        unsigned triggerID_old = trgnr_actual[i];
        // check for FEI4 problem, do not increase trigger to not loose synchronization
        // many correct triggers seem to be sent after the high trigger number storm (checked for DEPFET 2019 run 002141)
        unsigned cur_trigID = eudaq::PluginManager::GetTriggerID(*subevents[i].front());
        if (typeID[i] != tlu_id && cur_trigID > 32768) {
            std::cout << eudaq::Event::id2str(typeID[i]) << " " << subtype[i] << " " << cur_trigID << std::endl;
            // trgnr_actual[i]++;
            trgnr_out_of_order[i] = true;

        // catch NI error and increase NI total trigger number by one instead
        } else if (subtype[i] == "NI" &&
                   std::abs((long) (cur_trigID + trgnr_tlu_offset[i]) - (long) triggerID_old) > allowed_trg_jump_NI &&
                   !((trgnr_actual[i] + 1) % 32768 == 0 && !trgnr_out_of_order[i]) ) {
            std::cout << "NI problem " << cur_trigID + trgnr_tlu_offset[i] << " " << triggerID_old << std::endl;
            trgnr_actual[i] = (trgnr_actual[i] + 1) % 32768 + trgnr_tlu_offset[i];
            trgnr_out_of_order[i] = true;
            check_tlu_overflow(i, triggerID_old);
        } else {
            trgnr_actual[i] = trgnr_tlu_offset[i] + cur_trigID;
            // check for overflow                   // except if NI comes back from trigger problem in this event
            // if(!(trgnr_out_of_order[i] && subtype[i] == "NI"))
            check_tlu_overflow(i, triggerID_old);
            // reset trg_nr_problem flag
            if (trgnr_out_of_order[i])
                trgnr_out_of_order[i] = false;
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
        for (const auto &i : taglist) {
            events.back().SetTag(i, read_event->GetTag(i));
        }
        // read in sub events
        for (size_t i = 0; i < n_sub; i++) {
            // create shared pointer references to subevents in queues
            subevents[i].emplace(std::shared_ptr<eudaq::Event>(read_event->GetEventPtr(i)));
            // if queue was empty, determine trigger number of new current front element
            if (subevents[i].size() == 1) {
                update_trigger_id(i);
            }
        }
    }

    // checks for FEI4 high trigger number problem
    bool check_high_trg() {
        bool any_trgnr_out_of_order = false;
        for (size_t i = 0; i < n_sub; i++)
            any_trgnr_out_of_order |= trgnr_out_of_order[i];
        if (any_trgnr_out_of_order) {
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
            if ((typeID[i] != tlu_id) && (trgnr_actual[i] > trgnr_actual_global) && (trgnr_actual[i] < trgnr_actual_global + 10000) &&
                !trgnr_out_of_order[i]) {
                trgnr_actual_global = trgnr_actual[i];
                // mismatched_det = i;
            }
        }
    }

    // checks event for mismatch
    bool check_event() {
        // check trigger IDs for mismatch
        for (size_t i = 0; i < n_sub; i++) {
            if (typeID[i] != tlu_id && (trgnr_actual[i] != events.front().GetEventNumber() + 1 ||
                                         trgnr_out_of_order[i])) {  // && !subevents[i].front()->IsEORE()) {
                std::cout << i << " " << eudaq::Event::id2str(typeID[i]) << ":" << subtype[i] << " "
                          << trgnr_actual[i] << std::endl;
                // hack because NI sends invalid trigger number in EORE
                if (!subevents[i].front()->IsEORE())
                    return true;
            }
        }
        return false;
    }

    // discard global events and subevents of lower trigger IDs
    void discard_mismatch() {
        std::cout << "max_trigger_ID " << trgnr_actual_global << std::endl;
        std::cout << events.front() << std::endl;
        for (size_t i = 0; i < n_sub; i++)
            std::cout << *subevents[i].front() << std::endl;
        // check global event
        if ((events.front().GetEventNumber() + 1) != trgnr_actual_global) events.pop();
        // check sub events
        for (size_t i = 0; i < n_sub; i++) {
            // check tlu event
            if (typeID[i] == tlu_id) {
                if ((subevents[i].front()->GetEventNumber() + 1) != trgnr_actual_global) {
                    std::cout << "delete TLU subevent " << subevents[i].front()->GetEventNumber() << std::endl;
                    subevents[i].pop();
                }
            }
                // check other sub events
            else if (trgnr_actual[i] != trgnr_actual_global || trgnr_out_of_order[i]) {
                std::cout << "delete " << i << " " << eudaq::Event::id2str(typeID[i]) << ":" << subtype[i]
                          << " " << trgnr_actual[i] << std::endl;
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
};


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
        TFile rootfile((output_file.substr(0, output_file.length() - 3) + "root").c_str(), "RECREATE");
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

//        unsigned subevent_trigger[setup.n_sub];
        // for (size_t i = 0; i < n_sub; i++) event_tree.Branch("subevent_trigger", &subevent_trigger[i], "subevent_trigger/i");
        // event_tree.Branch("subevent_trigger", subevent_trigger, ("subevent_trigger["+std::to_string(setup.n_sub)+"]/i").c_str());


        // mismatch variables
        bool got_mismatch;

        while (reader.NextEvent()) { //  && counter // < 131100) { //  && counter < 114100) { // && counter < 37000) {
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
        while (setup.all_queues_filled()) {
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

        std::cout << "Remaining global events: " << setup.events.size() << std::endl;
        for (size_t i = 0; i < setup.n_sub; i++) {
            std::cout << "Remaining subevents: " << std::endl;
            std::cout << eudaq::Event::id2str(setup.typeID[i]) << " " << setup.subtype[i] << " "
                      << setup.subevents[i].size()
                      << std::endl;
        }
        event_tree.Write();
        rootfile.Close();

    } catch (...) {
        return op.HandleMainException();
    }

    return 0;
}
