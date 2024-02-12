#include <vector>
#include <map>
#include <tuple>
#include <iostream>
#include <cmath>
#include <queue>
#include <functional>
#include <algorithm>
#include <string>
#include <sstream>

// Define a custom type for the flow data
using SubnetData = std::tuple<int, int, int, int>; // qos, id, hop, flowId

using TimeData = std::tuple<int, int, int, int>; // tsai, tdi, tstcd, tci

using FlowData = std::tuple<int, int, int, int, int, int>; // hop, id, qos, k, tsai, flowId

// Define the maps to hold vectors of flows and times
std::map<std::string, std::vector<FlowData>> egressFlowset;
std::map<std::string, std::vector<TimeData>> egressTimeset;
std::map<std::string, std::vector<FlowData>> Flowset;
std::map<std::string, std::vector<TimeData>> Timeset;

using SubNetworkType = std::map<std::string, std::vector<SubnetData>>;

std::map<int, std::vector<std::vector<int>>> allocationTS;
std::map<int, std::vector<std::vector<int>>> talkerResults;
std::map<int, std::vector<std::string>> GCL;

// Assuming the existence of these types based on the Python code
using FlowIdRecord = std::map<int, int>;
using ShortestPathType = std::string;
using SrcDstType = std::vector<int>;
using QosType = std::vector<int>;
using ReverseLookupType = std::map<int, int>;
using FlowSrcDstType = std::map<int, std::vector<int>>;

SubNetworkType subNetwork;
FlowSrcDstType flowSrcDst;
	
// Assuming userDefined is a global or externally provide4d variable
bool userDefined = true; // You can change this to false to test the other condition



SubNetworkType TopologyPartition(const ShortestPathType& ShortestPath, const SrcDstType& SrcDst, int Pos, const QosType& Qos, const ReverseLookupType& reverse_lookup, const FlowIdRecord& flowIdRecord) {
    
    int lenLista = ShortestPath.size();
    int flowId = flowIdRecord.at(Pos);
    int srcEs = ShortestPath[0];
    int dstEs = ShortestPath[lenLista - 1];

    flowSrcDst[flowId].push_back(srcEs);
    flowSrcDst[flowId].push_back(dstEs);

    if (lenLista == 2) {
        subNetwork[ShortestPath].push_back(std::make_tuple(Qos[Pos], Pos + 1, 1, flowId));
        return subNetwork;
    }

    int startIndex = 0;
    int lenDictb = reverse_lookup.size();
    for (int i = 0; i < ShortestPath.size(); ++i) {
        int x = ShortestPath[i];
        if (i == 0) {
            continue;
        }
        if (reverse_lookup.find(x) != reverse_lookup.end()) {
            int qos;
            if (lenLista == i + 1 - startIndex) {
                qos = Qos[Pos];
            } else {
                qos = Qos[Pos] * (i + 1 - startIndex) / (lenLista + 1);
            }
            ShortestPathType pathSegment(ShortestPath.begin() + startIndex, ShortestPath.begin() + i + 1);
            subNetwork[pathSegment].push_back(std::make_tuple(qos, Pos + 1, startIndex + 1, flowId));
            startIndex = i;
            --lenDictb;
        }
        if (lenDictb == 0 && i < lenLista) {
            int qos;
            if (lenLista == lenLista - i) {
                qos = Qos[Pos];
            } else {
                qos = Qos[Pos] * (lenLista - i) / (lenLista + 1);
            }
            ShortestPathType pathSegment(ShortestPath.begin() + i, ShortestPath.end());
            subNetwork[pathSegment].push_back(std::make_tuple(qos, Pos + 1, i + 1, flowId));
            break;
        }
    }
    return subNetwork;
}


std::pair<int, int> DelayCalculation(int HopCnt) {
    const int FrameLength = 1500; // This is constrained by 's algorithm
    const int LinkSpeed = 1000;
    const int LinkLength = 100;
    if (userDefined) {
        int d_trans = (FrameLength + 30) * 8000 / LinkSpeed; // Unit is nano-second
        int IFG = 96000 / LinkSpeed; // Unit is nano-second
        int LenTSinTCI = d_trans + IFG;
        // Add a function to calculate pre_e2e
        int d_proc = 8000; // Unit is nano-second
        int d_prop = LinkLength * 10 / 2; // Unit is nano-second
        int arriveTSatCurrentSW = HopCnt * (d_proc + d_prop + d_trans);
        return std::make_pair(LenTSinTCI, arriveTSatCurrentSW);
    } else {
        return std::make_pair(12336, HopCnt * 20740);
    }
}

void assignToEgress(const std::map<std::string, std::vector<FlowData>>& Flowset,
const std::map<std::string, std::vector<TimeData>>& Timeset) {
    
    for (const auto& pair : Flowset) {
        const std::string& sn = pair.first;
        const std::vector<FlowData>& flow = pair.second;

        size_t lensn = sn.length();
        if (lensn == 2) {
            for (const auto& f : flow) {
                egressFlowset[sn].push_back(f);
            }
            egressTimeset[sn].push_back(Timeset.at(sn)[0]);
        } else {
            size_t numEgress = lensn - 1;
            std::string egress;
            for (size_t cnt = 0; cnt < lensn - 1; ++cnt) {
                egress = {sn[cnt], sn[cnt + 1]};
                for (const auto& f : flow) {
                    // Create a modified flow data with incremented hop
                    FlowData modifiedFlow = f;
                    std::get<0>(modifiedFlow) += cnt;
                    egressFlowset[egress].push_back(modifiedFlow);
                }
                egressTimeset[egress].push_back(Timeset.at(sn)[0]);
            }
        }
    }
}

void flowVarsForEgress(const std::map<std::string, std::vector<FlowData>>& Flowset) {
    
    for (const auto& pair : Flowset) {
        const std::string& sn = pair.first;
        const std::vector<FlowData>& flow = pair.second;

        size_t lensn = sn.length();
        if (lensn == 2) {
            for (const auto& f : flow) {
                egressFlowset[sn].push_back(f);
            }
        } else {
            size_t numEgress = lensn - 1;
            std::string egress;
            for (size_t cnt = 0; cnt < lensn - 1; ++cnt) {
                egress = {sn[cnt], sn[cnt + 1]};
                for (const auto& f : flow) {
                    egressFlowset[egress].push_back(f);
                }
            }
        }
    }
}

std::vector<int> step1andstep2_CalculateRelatedParameters(const SubNetworkType& subnetwork) {
    int numFlow = 0;
    for (const auto& entry : subnetwork) {
        const std::string& sn = entry.first;
        auto flow = entry.second; // Copy the flow set

        // Convert flow list to flow heap
        std::priority_queue<SubnetData, std::vector<SubnetData>, std::greater<SubnetData>> flowHeap(flow.begin(), flow.end());

        // Assuming flowHeap is not empty and has at least one element
        int tdi = std::get<0>(flowHeap.top()); // Get the smallest qos to be TDI

        int lenFlow = flow.size();
        double atstcd = 0.0; // The average number of time slot allocated in TCI

        while (!flowHeap.empty()) {
            SubnetData flowData = flowHeap.top();
            flowHeap.pop();
            int qos = std::get<0>(flowData);
            int id = std::get<1>(flowData);
            int hop = std::get<2>(flowData);
            int flowId = std::get<3>(flowData);

            if (numFlow <= id) {
                numFlow = id;
            }

            int x = qos / tdi; // Floor division
            int k = static_cast<int>(std::pow(2, std::floor(std::log2(x))));
            int tsai = k * tdi;
            atstcd += static_cast<double>(tdi) / tsai;

            
            Flowset[sn].push_back({hop, id, qos, k, tsai, flowId});
            
            if (!flow.empty()) { // All flows are traversed
	            int tstcd = std::ceil(atstcd);
	            int tci = tstcd * 12336;
	            Timeset[sn].push_back({tsai, tdi, tstcd, tci}); // At this time tsai is the max(tsai), thus is equal to the hyper period
	        }
        }

        
    }

    // Flowset contains all flow information, Timeset contains all time-related information for an egress
    assignToEgress(Flowset, Timeset); // Sub-network might contain many switches, this function breaks subNetwork to two adjacency switch set

    std::vector<int> FMTItalker(numFlow, 1);
    return FMTItalker;
}

void clearVars() {
    // Assuming flowSrcDst, talkerResults, TStoflow, GCL are defined elsewhere
    flowSrcDst.clear();
    talkerResults.clear();
    allocationTS.clear();
    egressFlowset.clear();
    egressTimeset.clear();
    Flowset.clear();
    Timeset.clear();
    subNetwork.clear();
    //TStoflow.clear();
    GCL.clear();
}


void GenerateGCLs(int ArrivedTimeInstance, int tdi, int numTDI, int tstcd, int egress, const std::string& gclname) {
    const int GB_interval = 12240;
    const std::string GateState_TSN = "10000000";
    const std::string GateState_BE = "01111111";
    const std::string GateState_GB = "00000000";

    int GCL_len = numTDI * 3;
    int TSN_interval = tstcd * 12336;
    int BE_interval = tdi - TSN_interval - GB_interval;

    int first_interval = ArrivedTimeInstance - GB_interval;
    int last_interval = BE_interval - first_interval;

    std::ostringstream gclStream;
    gclStream << "hyperPeriod    " << tdi * numTDI << "\n"
              << "numOfEntries    " << GCL_len << "\n"
              << GateState_BE << "    " << first_interval << "\n";

    for (int i = 1; i < numTDI; ++i) {
        gclStream << GateState_GB << "    " << GB_interval << "\n"
                  << GateState_TSN << "    " << TSN_interval << "\n"
                  << GateState_BE << "    " << BE_interval << "\n";
    }

    gclStream << GateState_GB << "    " << GB_interval << "\n"
              << GateState_TSN << "    " << TSN_interval << "\n"
              << GateState_BE << "    " << last_interval << "\n";

    // Add the generated GCL string to the GCL map
    GCL[egress].push_back(gclStream.str());
}

std::vector<int> convertStrToVector(const std::string& strdata, int flag) {
    std::vector<int> element;
    std::istringstream iss(strdata);
    std::string token;

    if (flag == 1) {
        while (getline(iss, token, ',')) {
            if (!token.empty() && token[0] != '[' && token[token.size() - 1] != ']') {
                element.push_back(std::stoi(token));
            }
        }
    } else if (flag == 2) {
        std::vector<std::vector<int>> subelement;
        while (getline(iss, token, ']')) {
            if (!token.empty()) {
                std::istringstream inner_iss(token.substr(1)); // Skip the '['
                std::vector<int> inner_element;
                std::string inner_token;
                while (getline(inner_iss, inner_token, ',')) {
                    if (!inner_token.empty()) {
                        inner_element.push_back(std::stoi(inner_token));
                    }
                }
                subelement.push_back(inner_element);
            }
        }
        // Assuming you want to flatten the 2D vector into a 1D vector
        for (const auto& vec : subelement) {
            element.insert(element.end(), vec.begin(), vec.end());
        }
    }
    return element;
}


int main() {
	
/*** for test TopologyPartition ***/
    // Example data for the function call
    ShortestPathType ShortestPath = {1, 2, 3, 4};
    SrcDstType SrcDst = {1, 4}; // Example source and destination
    int Pos = 0; // Example position
    QosType Qos = {10}; // Example QoS values
    ReverseLookupType reverse_lookup = {{2, 1}, {3, 2}}; // Example reverse lookup
    FlowIdRecord flowIdRecord = {{0, 1}}; // Example flow ID record

    // Call the TopologyPartition function with the example data
    subNetwork = TopologyPartition(ShortestPath, SrcDst, Pos, Qos, reverse_lookup, flowIdRecord);

    // Output the results (for demonstration purposes)
    for (const auto& pair : subNetwork) {
        std::cout << "Path: ";
        for (int node : pair.first) {
            std::cout << node << " ";
        }
        std::cout << "=> Values: ";
        for (const auto& tuple : pair.second) {
            std::cout << "(" << std::get<0>(tuple) << ", " << std::get<1>(tuple) << ", " << std::get<2>(tuple) << ", " << std::get<3>(tuple) << ") ";
        }
        std::cout << std::endl;
    }

/*** for test DelayCalcuation ***/
    int HopCnt = 5; // Example Hop Count
    std::pair<int, int> result = DelayCalculation(HopCnt);

    std::cout << "LenTSinTCI: " << result.first << ", arriveTSatCurrentSW: " << result.second << std::endl;


    return 0;
}
