//  main.cpp
//
//  Original version created by Grayson Rich 2014.
//  Version 2 created by Sam Hedges 2019.
//  This version created by Charlie Prior 2025.
//
//
//  modeled loosely after code by K. Kazkaz (LLNL)
//  original decoder by Kazkaz was used for 3320-generated data
//
//  this converter is meant to translate an NGM-binary file into a ROOT version
//  no analysis is performed, no data is removed except for NGM headers

//#include <stdio.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <cstdint>

#include <argparse/argparse.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TTreeIndex.h"

using s_clock = std::chrono::steady_clock;
namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
    fs::path inputFilename;
    fs::path outfileName;
    bool debug = false;
    bool clobber = false;
    bool quiet = false;
    int numCards = 1;
    int channelsPerCard = 16;

    argparse::ArgumentParser program("ngm2root");
    program.add_argument("input")
            .help("input file name")
            .store_into(inputFilename);

    program.add_argument("-o", "--output")
            .help("output file name or directory");

    program.add_argument("-n", "--numCards")
            .help("number of cards")
            .default_value(numCards)
            .scan<'i', int>()
            .store_into(numCards);

    program.add_argument("-c", "--channels")
            .help("channels per card")
            .default_value(channelsPerCard)
            .scan<'i', int>()
            .store_into(channelsPerCard);

    program.add_argument("--debug")
            .help("enable debug mode")
            .default_value(false)
            .implicit_value(true)
            .store_into(debug);

    program.add_argument("-q, ""--quiet")
            .help("enable quiet mode")
            .default_value(false)
            .implicit_value(true)
            .store_into(quiet);

    program.add_argument("-c", "--clobber")
            .help("overwrite output file if it exists")
            .default_value(false)
            .implicit_value(true)
            .store_into(clobber);

    try {
        program.parse_args(argc, argv);
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        std::cerr << program;
        return 1;
    }

    // check to see if the input file ends in '.bin' extension
    // if it doesn't, exit the program - user probably entered something wrong
    if (inputFilename.extension() != ".bin") {
        std::cerr << "Input file should end in .bin... did you do something wrong?" << std::endl;
        std::cerr << program;
        return 1;
    }

    if (auto out = program.present("-o")) {
        fs::path path(*out);
        if (fs::is_symlink(path)) {
            path = fs::read_symlink(path);
        }

        if (fs::is_directory(path)) {
            if (!fs::exists(path)) {
                std::cerr << "Output directory does not exist: " << path << std::endl;
                return 1;
            }
            // if the output path is a directory, generate the output file name
            // by appending the input file name to the output path
            outfileName = path / inputFilename.filename().replace_extension(".root");
        } else {
            // check parent directory exists
            if (!fs::exists(path.parent_path())) {
                std::cerr << "Parent directory of output file does not exist: " << path.parent_path() << std::endl;
                return 1;
            }
            outfileName = path;
        }
    } else {
        // if the output path is not specified, generate the output file name
        // by replacing the input file extension with .root
        outfileName = inputFilename.replace_extension(".root");
    }

    if (!clobber && fs::exists(outfileName)) {
        std::cerr << "Output file already exists: " << outfileName << std::endl;
        std::cerr << "Use -c or --clobber to overwrite." << std::endl;
        return 1;
    }

    if (!quiet) {
        std::cout << "Writing file to " << outfileName << std::endl;
    }

    std::string baseName = inputFilename.filename().replace_extension().string();

    //	Create the output file and TTree
    uint64_t timestamp;
    uint16_t peakHighIndex;
    uint16_t peakHighValue;
    uint16_t channelID;
    uint8_t formatBits;
    uint32_t accumulatorSum[8];
    uint8_t informationBits;
    uint32_t mawMaximumValue;
    uint32_t mawValueAfterTrigger;
    uint32_t mawValueBeforeTrigger;
    uint32_t mawTestData;
    uint32_t startEnergyValue;
    uint32_t maxEnergyValue;
    bool mawTestFlag;
    bool pileupFlag;
    uint32_t nSamples;
    uint16_t waveform[65536]; // max size of waveform (64K)

    uint32_t readBuffer[4096];
    char *bufferPointer = reinterpret_cast<char*>(readBuffer);
    uint32_t tmpWord;

    std::unique_ptr<TFile> outFile(TFile::Open(outfileName.c_str(), "RECREATE", baseName.c_str()));
    auto unsortedTree = new TTree("unsortedTree", baseName.c_str());
    unsortedTree->SetDirectory(nullptr);

    unsortedTree->Branch("ch", &channelID, "ch/s");
    unsortedTree->Branch("ts", &timestamp, "ts/l");
    unsortedTree->Branch("phi", &peakHighIndex, "phi/s");
    unsortedTree->Branch("phv", &peakHighValue, "phv/s");
    unsortedTree->Branch("fb", &formatBits, "fb/b");
    unsortedTree->Branch("acc", accumulatorSum, "acc[8]/i");
    unsortedTree->Branch("info", &informationBits, "info/b");
    unsortedTree->Branch("mawMax", &mawMaximumValue, "mawMax/i");
    unsortedTree->Branch("mawVat", &mawValueAfterTrigger, "mawVat/i");
    unsortedTree->Branch("mawVbt", &mawValueBeforeTrigger, "mawVbt/i");
    unsortedTree->Branch("mawTd", &mawTestData, "mawTd/i");
    unsortedTree->Branch("startEv", &startEnergyValue, "startEv/i");
    unsortedTree->Branch("maxEv", &maxEnergyValue, "maxEv/i");
    unsortedTree->Branch("mawTest", &mawTestFlag, "mawTest/O");
    unsortedTree->Branch("pu", &pileupFlag, "pu/O");
    unsortedTree->Branch("ns", &nSamples, "ns/i");
    unsortedTree->Branch("wf", waveform, "wf[ns]/s");

    auto finalTree = unsortedTree->CloneTree(0);
    finalTree->SetName("t");
    finalTree->SetDirectory(outFile.get());

    uint64_t index = 0;
    const auto processStartingTime = s_clock::now(); // overall timer start
    uint32_t packetWords; // number of words in packet for channel
    bool incompleteDumpFlag = false;
    uint64_t spillNumber = 0;

    // open the input file
    std::ifstream inFile(inputFilename, std::ios::binary);

    //Get size of file, from:
    //https://stackoverflow.com/questions/2409504/using-c-filestreams-fstream-how-can-you-determine-the-size-of-a-file
    inFile.ignore(std::numeric_limits<std::streamsize>::max());
    std::streamsize bitsInFile = inFile.gcount();
    inFile.clear();
    inFile.seekg(0, std::ios_base::beg);

    if (!quiet) { std::cout << "Processing " << inputFilename << "..." << std::endl; }

    // seek past the header at the start of the binary file
    // 100 words = 400 bytes
    inFile.seekg(400);
    bitsInFile -= 400;

    if (!quiet) { std::cout << "Header seeked through.." << std::endl; }

    while (inFile.good()) {
        // we're at the start of a spill
        // read 10 word spill header
        inFile.read(bufferPointer, 40);
        bitsInFile -= 40;

        if (debug) {
            // if we're debugging, print out the first 10 words of the first spill
            for (int i = 0; i < 10; i++) {
                memcpy(&tmpWord, &readBuffer[i], 4);
                if (!quiet) {
                    std::cout << "Word " << i << " of " << spillNumber << " spill: \t "
                            << std::hex << tmpWord << std::dec << std::endl;
                }
            }
        }

        // check for EOF
        memcpy(&tmpWord, &readBuffer[0], 4);
        if (tmpWord == 0x0E0F0E0F) {
            break;
        }

        for (int cardNumber = 0; cardNumber < numCards; cardNumber++) {
            // skip the packet header for the spill
            // this is two words
            inFile.read(bufferPointer, 8);
            bitsInFile -= 8;

            // now that we're inside the spill, we have to parse data for each channel
            for (int channelNumber = 0; channelNumber < channelsPerCard; channelNumber++) {
                // read in the packet header for the channel
                inFile.read(bufferPointer, 32);
                bitsInFile -= 32;

                // from the channel info packet header, we can determine the size of the packet data
                memcpy(&tmpWord, &readBuffer[7], 4);
                packetWords = tmpWord;
                bitsInFile -= (packetWords * 4);


                if (debug && !quiet) {
                    std::cout << "Bits in file before digitizer dump: " << bitsInFile + packetWords * 4 << std::endl;
                    std::cout << "Number of words in packet for channel " << channelNumber << " in spill "
                            << spillNumber << ": \t " << packetWords << std::endl;
                    std::cout << "Bits in file after digitizer dump: " << bitsInFile << std::endl;
                }

                //Check if more words are expected than are in file:
                if (bitsInFile < 40 && !quiet) {
                    incompleteDumpFlag = true;
                    std::cout << "Found incomplete dump, remaining data not added to TTree" << std::endl;
                    break;
                }

                while (packetWords > 0) {
                    if (index % 100000 == 0 && !quiet) {
                        const auto endTime = s_clock::now();
                        const std::chrono::duration<double> elapsedTime = endTime - processStartingTime;
                        const auto elapsedTimeSeconds = std::chrono::duration_cast<std::chrono::seconds>(elapsedTime);
                        std::cout << "Processed " << index << " events in "
                                << std::setprecision(0) << elapsedTimeSeconds.count() << std::setprecision(6)
                                << " seconds" << std::endl;
                    }

                    // first two words of an event are there no matter what the format bits are set to
                    inFile.read(bufferPointer, 8); packetWords -= 2; memcpy(&tmpWord, &readBuffer[0], 4);
                    formatBits = static_cast<uint8_t>(tmpWord & 0xf);
                    channelID = static_cast<uint16_t>((tmpWord & 0xfff0) >> 4);
                    timestamp = (static_cast<uint64_t>(tmpWord & 0xffff0000)) << 16;

                    memcpy(&tmpWord, &readBuffer[1], 4);
                    timestamp = timestamp | tmpWord;

                    if (debug && !quiet) {
                        // print out the first two words of the event
                        memcpy(&tmpWord, &readBuffer[0], 4);
                        std::cout << "First two words of event:\t " << std::hex << tmpWord << "\t";
                        memcpy(&tmpWord, &readBuffer[1], 4);
                        std::cout << tmpWord << std::dec << std::endl;
                    }

                    if (debug && !quiet) {
                        // print out the determined format bits
                        std::cout << "Format bits: " << std::bitset<4>(formatBits) << std::endl;
                    }

                    if ((formatBits & 0x1) != 0) {
                        if (debug && !quiet) {
                            std::cout << "Reading words for format bit 0" << std::endl;
                        }

                        inFile.read(bufferPointer, 7 * 4); packetWords -= 7; memcpy(&tmpWord, &readBuffer[0], 4);
                        peakHighValue = tmpWord & 0xffff;
                        peakHighIndex = (tmpWord & 0xffff0000) >> 16;

                        if (debug && !quiet) {
                            std::cout << "Peak index: " << peakHighIndex
                                    << " peak value: " << peakHighValue << std::endl;
                        }

                        memcpy(&tmpWord, &readBuffer[1], 4);
                        informationBits = (tmpWord & 0xff000000) >> 24;
                        accumulatorSum[0] = (tmpWord & 0xffffff);

                        memcpy(&accumulatorSum[1], &readBuffer[2], 4);
                        memcpy(&accumulatorSum[2], &readBuffer[3], 4);
                        memcpy(&accumulatorSum[3], &readBuffer[4], 4);
                        memcpy(&accumulatorSum[4], &readBuffer[5], 4);
                        memcpy(&accumulatorSum[5], &readBuffer[6], 4);
                    } else {
                        peakHighIndex = 0;
                        peakHighValue = 0;
                        informationBits = 0;
                        accumulatorSum[0] = 0;
                        accumulatorSum[1] = 0;
                        accumulatorSum[2] = 0;
                        accumulatorSum[3] = 0;
                        accumulatorSum[4] = 0;
                        accumulatorSum[5] = 0;
                    }
                    if ((formatBits & 0x2) != 0) {
                        if (debug && !quiet) {
                            std::cout << "Reading words for format bit 1" << std::endl;
                        }
                        inFile.read(bufferPointer, 2 * 4); packetWords -= 2;
                        memcpy(&accumulatorSum[6], &readBuffer[0], 4);
                        memcpy(&accumulatorSum[7], &readBuffer[1], 4);
                    } else {
                        // populate accumulators with 0 if they're not defined
                        accumulatorSum[6] = 0;
                        accumulatorSum[7] = 0;
                    }
                    if ((formatBits & 0x4) != 0) {
                        if (debug && !quiet) {
                            std::cout << "Reading words for format bit 2" << std::endl;
                        }
                        inFile.read(bufferPointer, 3 * 4); packetWords -= 3;
                        memcpy(&mawMaximumValue, &readBuffer[0], 4);
                        memcpy(&mawValueAfterTrigger, &readBuffer[1], 4);
                        memcpy(&mawValueBeforeTrigger, &readBuffer[2], 4);
                    } else {
                        mawMaximumValue = 0;
                        mawValueAfterTrigger = 0;
                        mawValueBeforeTrigger = 0;
                    }
                    if ((formatBits & 0x8) != 0) {
                        if (debug && !quiet) {
                            std::cout << "Reading words for format bit 3" << std::endl;
                        }
                        inFile.read(bufferPointer, 2 * 4); packetWords -= 2;
                        memcpy(&startEnergyValue, &readBuffer[0], 4);
                        memcpy(&maxEnergyValue, &readBuffer[1], 4);
                    } else {
                        startEnergyValue = 0;
                        maxEnergyValue = 0;
                    }

                    // the next word will determine the number of sample words we read
                    inFile.read(reinterpret_cast<char*>(&tmpWord), 4); packetWords -= 1;
                    nSamples = 2 * (tmpWord & 0x3ffffff);

                    if (debug && !quiet) {
                        std::cout << "Determined there are " << nSamples << " sample words" << std::endl;
                    }

                    pileupFlag = (tmpWord & 0x4000000) >> 26;
                    mawTestFlag = (tmpWord & 0x8000000) >> 27;

                    if (mawTestFlag) {
                        std::cerr << "MAW test flag is set, but this is not supported" << std::endl;
                        return 1;
                    }

                    for (int i = 0; i < (nSamples / 2); i++) {
                        inFile.read(reinterpret_cast<char*>(&tmpWord), 4); packetWords -= 1;
                        waveform[i * 2] = tmpWord & 0xffff;
                        waveform[i * 2 + 1] = (tmpWord & 0xffff0000) >> 16;
                    }

                    mawTestData = 0; // Don't support reading MAW test data

                    unsortedTree->Fill();
                    index++;
                }
            }
        }


        if (!incompleteDumpFlag) {
            // now we do the time ordering of the spill
            const auto startTime = s_clock::now();

            if (debug && !quiet) {
                std::cout << "Creating TTree index for spill " << spillNumber << std::endl;
            }

            unsortedTree->BuildIndex("ts");
            auto treeIndex = reinterpret_cast<TTreeIndex *>(unsortedTree->GetTreeIndex());

            if (debug && !quiet) {
                const auto endTime = s_clock::now();
                const std::chrono::duration<double> elapsedTime = endTime - startTime;;
                const auto elapsedTimeSeconds = std::chrono::duration_cast<std::chrono::seconds>(elapsedTime);
                std::cout << "Index created for spill " << spillNumber << " in "
                        << std::setprecision(0) << elapsedTimeSeconds.count() << std::setprecision(6)
                        << " seconds" << std::endl;
            }

            if (debug && !quiet) {
                std::cout << "Tree index contains " << treeIndex->GetN() << " entries" << std::endl;
            }

            for (int i = 0; i < treeIndex->GetN(); i++) {
                unsortedTree->GetEntry(treeIndex->GetIndex()[i]);

                if (debug && !quiet) {
                    std::cout << i << "-th event recalled is index " << treeIndex->GetIndex()[i] << std::endl;
                }

                finalTree->Fill();
            }


            spillNumber++;

            treeIndex->Delete();
            unsortedTree->Reset();
        }
    }

    unsortedTree->Delete();
    finalTree->Write("sis3316tree", TObject::kOverwrite);
    outFile->Close();
    inFile.close();

    if (!quiet) {
        const auto endTime = s_clock::now();
        const std::chrono::duration<double> elapsedTime = endTime - processStartingTime;;
        const auto elapsedTimeSeconds = std::chrono::duration_cast<std::chrono::seconds>(elapsedTime);
        std::cout << "Recorded " << index << " events in "
                << std::setprecision(0) << elapsedTimeSeconds.count() << std::setprecision(6)
                << " seconds" << std::endl;
    }

    return 0;
}
