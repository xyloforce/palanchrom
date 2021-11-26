#include <iostream>
#include <fstream>
#include <regex>

#include "fasta_tools.h"
#include "bed_tools.h"

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        throw std::domain_error("Unsufficient number of args : need source fasta and name of bed output");
    }

    std::cout << "Reading input..." << std::endl;
    fasta inputFile(argv[1], "read_line", false);

    std::cout << "Creating output..." << std::endl;
    bed outputFile(argv[2], false);

    std::cout << "Starting analysis..." << std::endl;
    int count(0);

    while (!inputFile.isEOF())
    {
        fasta_entry entry(inputFile.read_fasta_line());
        int startCG = 0;
        int stopCG = 0;

        std::string sequenceData(entry.getSequence());
        for (int i(0); i < sequenceData.size(); i++)
        {
            if (sequenceData[i] == 'C')
            {
                if (i - 2 > -1 && i + 1 < sequenceData.size())
                {
                    if (sequenceData.substr(i-2, 2) == "CG" && sequenceData.substr(i, 2) == "CG")
                    {
                        stopCG = i; // continuity of an int
                    } else if (sequenceData.substr(i, 2) == "CG") {
                        startCG = i;
                        stopCG = i; // new int
                    }
                }
            }
            else if (sequenceData[i] == 'G')
            {
                if (i - 1 > -1 && i + 2 < sequenceData.size())
                {
                    if (sequenceData.substr(i-1, 2) == "CG" && sequenceData.substr(i+1, 2) == "CG")
                    {
                        stopCG = i; // continuity of an int
                    }
                    else if (sequenceData.substr(i-1, 2) == "CG")
                    {
                        stopCG = i + 1;
                        outputFile.writeBedLine(bed_entry(entry.getChrom(), startCG, stopCG, ".", 0, '+'));
                    } // end of an int
                }
            }
            else
            {
            }
        }
        // char lastChar('\0');
        // for(int i(0); i<sequenceData.size(); i++) {
        //     if(lastChar != '\0') {
        //         if(toupper(lastChar) == 'C' && toupper(sequenceData[i]) == 'G') {
        //             cpg = cpg + lastChar + sequenceData[i];
        //             no_cpg += "--";
        //             lastChar = '\0';
        //         } else {
        //             cpg += "--";
        //             no_cpg = no_cpg + lastChar + sequenceData[i];
        //             lastChar = '\0';
        //         }
        //     } else if((i+1) > sequenceData.size()) {
        //         no_cpg += sequenceData[i];
        //     } else {
        //         lastChar = sequenceData[i];
        //     }
        count++;
        if (count % 10 == 0)
        {
            std::cout << "Treated " << count << " lines                      \r";
        }
    } // sequence done

    // write bed entry

    return 0;
}
