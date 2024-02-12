#include <vector>
#include <map>
#include <tuple>
#include <iostream>

// Define a custom type for the flow data
using FlowData = std::tuple<int, int, int, int, int, int>;
using TimeData = std::tuple<int, int, int, int>;

// Define the maps to hold vectors of flows and times
std::map<std::string, std::vector<FlowData>> egressFlowset;
std::map<std::string, std::vector<TimeData>> egressTimeset;
std::map<std::string, std::vector<FlowData>> Flowset;
std::map<std::string, std::vector<TimeData>> Timeset;

using SubNetworkType = std::map<std::vector<int>, std::vector<std::tuple<int, int, int, int>>>;

// Assuming the existence of these types based on the Python code
using FlowIdRecord = std::map<int, int>;
using ShortestPathType = std::vector<int>;
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
