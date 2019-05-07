//
// Created by paschen on 3/7/19.
//

#include <iostream>
#include <string>
#include "eudaq/OptionParser.hh"
#include "eudaq/FileReader.hh"
#include "eudaq/FileWriter.hh"
#include "eudaq/PluginManager.hh"


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

        // reader instance
        eudaq::FileReader reader(infile.Value());
        // pointer to constant detector event
        eudaq::DetectorEvent const *read_event;

        // writer instance
        eudaq::FileWriter *writer = eudaq::FileWriterFactory::Create("native");
        writer->SetFilePattern(output_file);
        writer->StartRun(reader.RunNumber());

        // display file name
        std::cout << "Reading file: " << reader.Filename() << std::endl;
        // get first event (BORE - Beginning Of Run Event) and copy to output file
        read_event = &reader.GetDetectorEvent();
        std::cout << "Writing to file: " << output_file << std::endl;
        writer->WriteEvent(*read_event);
        // initialize PluginManager with BORE
        eudaq::PluginManager::Initialize(*read_event);

        // variables used in the loop
        unsigned counter = 0;
        unsigned end = 1000000;

        unsigned runnumber;
        unsigned eventnumber;
        uint64_t timestamp;

        // deque for global events with global information
        std::deque<eudaq::DetectorEvent> events;
        unsigned n_sub = read_event->NumEvents();
        // deques for subevents
        std::deque<std::shared_ptr<eudaq::Event>> subevents[n_sub];
        unsigned long triggerIDs[n_sub];
        unsigned long triggerID_offsets[n_sub];
        for (size_t i = 0; i < n_sub; i++) triggerIDs[i] = triggerID_offsets[i] = 0;
        unsigned typeIDs[n_sub];
        std::string subtype[n_sub];
        for (size_t i = 0; i < n_sub; i++) {
            typeIDs[i] = read_event->GetEvent(i)->get_id();
            subtype[i] = read_event->GetEvent(i)->GetSubType();
        }
        unsigned tlu_id = eudaq::Event::str2id("_TLU");  // ID of TLU subevents

        // mismatch variables
        bool got_mismatch = false;
        bool high_trg_nr_problem[n_sub];
        for (size_t i = 0; i < n_sub; i++) high_trg_nr_problem[i] = false;
        unsigned long max_triggerID = 0;
        size_t mismatched_det = 0;
        unsigned long triggerID_old;

        while (reader.NextEvent()) {
            counter++;
            // get current event as detector event
            read_event = &reader.GetDetectorEvent();
            // get run number, event number and timestamp for copying
            runnumber = read_event->GetRunNumber();
            eventnumber = read_event->GetEventNumber();
            timestamp = read_event->GetTimestamp();

            // create global event in deque
            // FIXME: Replace by emplace_back to prevent temporary copy?
            events.push_back(eudaq::DetectorEvent(runnumber, eventnumber, timestamp));
            // copy flags
            events.back().SetFlags(read_event->GetFlags());
            // copy tags
            std::vector<std::string> taglist = read_event->GetTagList("");
            for (size_t i = 0; i < taglist.size(); i++) {
                events.back().SetTag(taglist[i], read_event->GetTag(taglist[i]));
            }
            // read in sub events
            for (size_t i = 0; i < n_sub; i++) {
                // if deque empty, determine new trigger ID
                if (subevents[i].empty()) {
                    triggerID_old = triggerIDs[i];
                    // prevent loss of trigger ID synchronisation for FEI4 problem
                    if (typeIDs[i] != tlu_id &&
                        eudaq::PluginManager::GetTriggerID(*read_event->GetEventPtr(i)) > 32768) {
                        std::cout << eudaq::Event::id2str(typeIDs[i]) << " " << subtype[i] << " "
                                  << eudaq::PluginManager::GetTriggerID(*read_event->GetEventPtr(i)) << std::endl;
                        triggerIDs[i]++;
                        high_trg_nr_problem[i] = true;
                    } else {
                        if (high_trg_nr_problem[i])
                            high_trg_nr_problem[i] = false;
                        triggerIDs[i] =
                                triggerID_offsets[i] + eudaq::PluginManager::GetTriggerID(*read_event->GetEventPtr(i));
                        // check for overflow  TODO: define detection limit
                        if (triggerIDs[i] < std::max((long) triggerID_offsets[i] + 1, (long) triggerID_old - 32760) &&
                            triggerID_old - triggerIDs[i] > 32000) {
                            std::cout << "overflow detected " << subtype[i] << " new: " << triggerIDs[i] << " old: "
                                      << triggerID_old << std::endl;
                            triggerID_offsets[i] += 32768;
                            triggerIDs[i] =
                                    triggerID_offsets[i] +
                                    eudaq::PluginManager::GetTriggerID(*read_event->GetEventPtr(i));
                        }
                    }
                }
                // create shared pointer references to subevents in deques
                // FIXME: Replace by emplace_back to prevent temporary copy?
                subevents[i].push_back(std::shared_ptr<eudaq::Event>(read_event->GetEventPtr(i)));
            }

            // if FEI4 trigger number problem, immediately discard entire event
            bool any_high_trg_nr_problem = false;
            for (size_t i = 0; i < n_sub; i++)
                any_high_trg_nr_problem |= high_trg_nr_problem[i];
            if (any_high_trg_nr_problem) {
                std::cout << events.front() << std::endl;
                events.pop_front();
                for (size_t i = 0; i < n_sub; i++) {
                    std::cout << *subevents[i].front() << std::endl;
                    subevents[i].pop_front();
                }
                std::cout << "High trigger number problem, skip event" << std::endl << std::endl;
                continue;
            }

            // check trigger IDs for mismatch
            bool new_mismatch = false;
            for (size_t i = 0; i < n_sub; i++) {
                if ((typeIDs[i] != tlu_id) && (triggerIDs[i] != (events.front().GetEventNumber() + 1))) {
                    std::cout << i << " " << eudaq::Event::id2str(typeIDs[i]) << ":" << subtype[i] << " "
                              << triggerIDs[i] << std::endl;
                    got_mismatch = true;
                    new_mismatch = true;
                }
                // TODO: define max trigger jump
                if ((typeIDs[i] != tlu_id) && (triggerIDs[i] > max_triggerID) &&
                    (triggerIDs[i] < max_triggerID + 10000)) {
                    max_triggerID = triggerIDs[i];
                    mismatched_det = i;
                }
            }
            if (!new_mismatch) {
                if (got_mismatch)
                    std::cout << "reset mismatch" << std::endl;
                got_mismatch = false;
            }

            // in case of mismatch discard global events and subevents of lower trigger IDs
            if (got_mismatch) {
                std::cout << "max_trigger_ID " << max_triggerID << std::endl;
                std::cout << events.front() << std::endl;
                for (size_t i = 0; i < n_sub; i++)
                    std::cout << *subevents[i].front() << std::endl;
                // return 0;
                // check global event
                if ((events.front().GetEventNumber() + 1) != max_triggerID) events.pop_front();
                // check sub events
                for (size_t i = 0; i < n_sub; i++) {
                    // check tlu event
                    if (typeIDs[i] == tlu_id) {
                        if ((subevents[i].front()->GetEventNumber() + 1) != max_triggerID) {
                            std::cout << "delete TLU subevent " << subevents[i].front()->GetEventNumber() << std::endl;
                            subevents[i].pop_front();
                        }
                    }
                        // check other sub events
                    else if (triggerIDs[i] != max_triggerID) {
                        std::cout << "delete " << i << " " << eudaq::Event::id2str(typeIDs[i]) << ":" << subtype[i]
                                  << " " << triggerIDs[i] << std::endl;
                        subevents[i].pop_front();
                        if (!subevents[i].empty()) {
                            triggerID_old = triggerIDs[i];
                            triggerIDs[i] =
                                    triggerID_offsets[i] + eudaq::PluginManager::GetTriggerID(*subevents[i].front());
                            // check for overflow  TODO: define detection limit
                            if (triggerIDs[i] <
                                std::max((long) triggerID_offsets[i] + 1, (long) triggerID_old - 32760) &&
                                triggerID_old - triggerIDs[i] > 32000) {
                                std::cout << "overflow detected " << subtype[i] << " new: " << triggerIDs[i] << " old: "
                                          << triggerID_old << std::endl;
                                triggerID_offsets[i] += 32768;
                                triggerIDs[i] =
                                        triggerID_offsets[i] +
                                        eudaq::PluginManager::GetTriggerID(*subevents[i].front());
                            }
                        }
                    }
                }
            }
                // no mismatch: write out event
            else {
                for (size_t i = 0; i < n_sub; i++) {
                    events.front().AddEvent(subevents[i].front());
                    // delete subevents[i]->front();
                    subevents[i].pop_front();
                    if (!subevents[i].empty()) {
                        triggerID_old = triggerIDs[i];
                        triggerIDs[i] =
                                triggerID_offsets[i] + eudaq::PluginManager::GetTriggerID(*subevents[i].front());
                        // check for overflow  TODO: define detection limit
                        if (triggerIDs[i] < std::max((long) triggerID_offsets[i] + 1, (long) triggerID_old - 32760) &&
                            triggerID_old - triggerIDs[i] > 32000) {
                            std::cout << "overflow detected " << subtype[i] << " new: " << triggerIDs[i] << " old: "
                                      << triggerID_old << std::endl;
                            triggerID_offsets[i] += 32768;
                            triggerIDs[i] =
                                    triggerID_offsets[i] + eudaq::PluginManager::GetTriggerID(*subevents[i].front());
                        }
                    }
                }
                // write out global event in deque and delete
                writer->WriteEvent(events.front());
                events.pop_front();
            }
        }

        // TODO: Write out possibly remaining events

        std::cout << "Remaining global events: " << events.size() << std::endl;
        for (unsigned i = 0; i < n_sub; i++) {
            std::cout << "Remaining subevents: " << std::endl;
            std::cout << eudaq::Event::id2str(typeIDs[i]) << " " << subtype[i] << " " << subevents[i].size()
                      << std::endl;
        }
    } catch (...) {
        return op.HandleMainException();
    }

    return 0;
}
