#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <string>
#include <regex>

void calepsilon(
    const char* input_list = "/mnt/e/git-repo/nuclear-structure/output/method2/NeEQMD/inputfile.list",
    const char* outputfile = "/mnt/e/git-repo/nuclear-structure/analysis/anaout/NeNe_EQMD_all_epsilon2.root"
)
{
    std::ifstream infile(input_list);
    if (!infile.is_open()) {
        std::cerr << "Error opening input list file." << std::endl;
        return;
    }

    TFile *file_out = TFile::Open(outputfile, "RECREATE");
    if (!file_out || file_out->IsZombie()) {
        std::cerr << "Error creating output file." << std::endl;
        return;
    }

    TTree *tree_out = new TTree("tree_out", "tree_out");

    // Loop over each line in the input list
    std::vector<double> probabilities = {0.56353591, 0.2854512, 0.00368324, 0.13075506, 0.01657459};
    std::vector<std::string> suffix_list;
    std::vector<double> probability_list;
    std::vector<double> epsilon2_xyz_list;
    std::vector<double> epsilon2_2_list;
    std::vector<double> epsilon2_4_list;

    std::string inputfile;
    while (std::getline(infile, inputfile)) 
    {
        // Extract EQMDij from the input filename
        std::smatch match;
        std::regex regex("EQMD(\\d{2})");
        std::string branch_suffix;
        int index_i, index_j;

        if (std::regex_search(inputfile, match, regex)) 
        {
            branch_suffix = match.str(1);
            suffix_list.push_back(branch_suffix);

            index_i = branch_suffix[0] - '0'; // 转换字符为整数
            index_j = branch_suffix[1] - '0';
            double probability = probabilities[index_i - 1] * probabilities[index_j - 1];
            probability_list.push_back(probability);
        } 
        else 
        {
            std::cerr << "Error: Could not extract EQMDij from file name: " << inputfile << std::endl;
            continue;
        }

        TFile *file = TFile::Open(inputfile.c_str(), "READ");
        if (!file || file->IsZombie()) 
        {
            std::cerr << "Error opening file: " << inputfile << std::endl;
            if (file) file->Close();
            continue;
        }

        TTree *tree_epsilon = (TTree*)file->Get("tree_epsilon");
        if (!tree_epsilon) 
        {
            std::cerr << "Error: tree_epsilon not found in file: " << inputfile << std::endl;
            file->Close();
            continue;
        }

        double epsilon2_P;
        tree_epsilon->SetBranchAddress("epsilon2_xyz_array", &epsilon2_P);

        double sum_epsilon2_xyz = 0.0;
        double sum_epsilon2_P2 = 0.0;
        double sum_epsilon2_P4 = 0.0;

        // Loop over the entries to calculate sums
        int nEntries = tree_epsilon->GetEntries();
        for (int i = 0; i < nEntries; ++i) 
        {
            tree_epsilon->GetEntry(i);
            sum_epsilon2_xyz += epsilon2_P;

            double epsilon2_P2 = epsilon2_P * epsilon2_P;
            sum_epsilon2_P2 += epsilon2_P2;
            sum_epsilon2_P4 += epsilon2_P2 * epsilon2_P2;
        }
        double avg_epsilon2_xyz = sum_epsilon2_xyz / nEntries;
        double avg_epsilon2_P2 = sum_epsilon2_P2 / nEntries;
        double avg_epsilon2_P4 = sum_epsilon2_P4 / nEntries;

        double epsilon2_2 = TMath::Sqrt(avg_epsilon2_P2);
        double epsilon2_4 = TMath::Power((2 * avg_epsilon2_P2 * avg_epsilon2_P2 - avg_epsilon2_P4), 0.25);

        std::cout << "avg_epsilon2_xyz: " << avg_epsilon2_xyz << std::endl;
        std::cout << "avg_epsilon2_P2: " << avg_epsilon2_P2 << std::endl;
        std::cout << "avg_epsilon2_P4: " << avg_epsilon2_P4 << std::endl;
        std::cout << "epsilon2_2: " << epsilon2_2 << std::endl;
        std::cout << "epsilon2_4: " << epsilon2_4 << std::endl;

        epsilon2_xyz_list.push_back(avg_epsilon2_xyz);
        epsilon2_2_list.push_back(epsilon2_2);
        epsilon2_4_list.push_back(epsilon2_4);
    }

    // 检查向量长度是否一致
    if (suffix_list.size()!=epsilon2_2_list.size()||suffix_list.size()!=epsilon2_4_list.size()||suffix_list.size()!=epsilon2_xyz_list.size()) 
    {
        std::cerr << "Error: Mismatch in vector lengths! "
                  << "suffix_list size: " << suffix_list.size() 
                  << ", epsilon2_list size: " << epsilon2_2_list.size() 
                  << ", epsilon4_list size: " << epsilon2_4_list.size() 
                  << std::endl;
        return; // 退出函数
    }

    // 根据 suffix 动态创建 branch，并将 vector 中的数据写入 branch
    for (size_t i = 0; i < suffix_list.size(); ++i) 
    {
        std::string branch_name_probability = "probability_list_EQMD" + suffix_list[i];
        std::string branch_name_xyz = "epsilon2_xyz_EQMD" + suffix_list[i];
        std::string branch_name_2 = "epsilon2_2_EQMD" + suffix_list[i];
        std::string branch_name_4 = "epsilon2_4_EQMD" + suffix_list[i];

        tree_out->Branch(branch_name_probability.c_str(), &probability_list[i], (branch_name_probability + "/D").c_str());
        tree_out->Branch(branch_name_xyz.c_str(), &epsilon2_xyz_list[i], (branch_name_xyz + "/D").c_str());
        tree_out->Branch(branch_name_2.c_str(), &epsilon2_2_list[i], (branch_name_2 + "/D").c_str());
        tree_out->Branch(branch_name_4.c_str(), &epsilon2_4_list[i], (branch_name_4 + "/D").c_str());
    }
    tree_out->Fill();
    
    file_out->cd();
    tree_out->Write();
    file_out->Close();

    infile.close();
}
