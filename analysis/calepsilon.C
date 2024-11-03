#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <iostream>

void calepsilon
(
    const char* inputfile = "/mnt/e/git-repo/nuclear-structure/output/method2/NeEQMD/NeNe_EQMD11_method2_100k.root",
    const char* outputfile= "/mnt/e/git-repo/nuclear-structure/analysis/anaout/NeNe_EQMD_test_epsilon2.root"
)
{
    TFile *file = TFile::Open(inputfile, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file." << std::endl;
        return;
    }

    TTree *tree_epsilon = (TTree*)file->Get("tree_epsilon");
    double epsilon2_P;
    tree_epsilon->SetBranchAddress("epsilon2_xyz_array", &epsilon2_P);

    double sum_epsilon2_P2 = 0.0;
    double sum_epsilon2_P4 = 0.0;
    int nEntries = tree_epsilon->GetEntries();

    // Loop over the entries
    for (int i = 0; i < nEntries; ++i) {
        tree_epsilon->GetEntry(i);
        double epsilon2_P2 = epsilon2_P * epsilon2_P;
        sum_epsilon2_P2 += epsilon2_P2;
        sum_epsilon2_P4 += epsilon2_P2 * epsilon2_P2;
    }

    std::cout << "sum_Epsilon_P{2} = " << sum_epsilon2_P2 << std::endl;
    std::cout << "sum_Epsilon_P{4} = " << sum_epsilon2_P4 << std::endl;

    // Calculate the average values
    double avg_epsilon2_P2 = sum_epsilon2_P2 / nEntries;
    double avg_epsilon2_P4 = sum_epsilon2_P4 / nEntries;

    std::cout << "avg_epsilon2_P2 = " << avg_epsilon2_P2 << std::endl;
    std::cout << "avg_epsilon2_P4 = " << avg_epsilon2_P4 << std::endl;

    // Calculate epsilon_n{2} and epsilon_n{4}
    double epsilon2_2 = TMath::Sqrt(avg_epsilon2_P2);
    double epsilon2_4 = TMath::Power((2 * avg_epsilon2_P2 * avg_epsilon2_P2 - avg_epsilon2_P4), 0.25);

    // Print results
    std::cout << "Epsilon_n{2} = " << epsilon2_2 << std::endl;
    std::cout << "Epsilon_n{4} = " << epsilon2_4 << std::endl;

    // Close the file
    file->Close();
}

