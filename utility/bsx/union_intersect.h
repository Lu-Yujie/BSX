// used to compute set intersection, faster in specific scenario
//   multi short sets, intersect with a long set
// not used so far

#ifndef MULTI_UNION_INTERSECT_H
#define MULTI_UNION_INTERSECT_H

#include <vector>
#include <set>
#include <map>
#include <queue>

using namespace std;
// Counting all the comparison times is too complicated, loop conditions also need to be counted
// Just count the time directly, and let's talk about time improvements

struct Compare {
    bool operator()(const pair<int, int>& lhs, const pair<int, int>& rhs) const {
        return lhs.first > rhs.first;
    }
};

// slowest
vector<set<int>>
union_concur_intersect_map(const vector<set<int>>& unions, const set<int>& intersect) {
    vector<set<int>> results;
    results.resize(unions.size(), {});
    map<int, vector<int>> topns;  // used as queue, keep min values of each unions are in queue

    // init iterators
    vector<set<int>::iterator> uni_iters;
    for (int i = 0; i < unions.size(); i++) {
        auto& array = unions[i];
        int first_ele = *(array.begin());
        topns[first_ele].emplace_back(i);
        // pair<int, vector<int>> first_ele;
        // first_ele.first = *(array.cbegin());
        // auto topn = topns.find(first_ele.first);
        // if (topn != topns.end()) {
        //     topn->second.emplace_back(i);
        // } else {
        //     first_ele.second.reserve(unions.size());
        //     first_ele.second.emplace_back(i);
        //     topns.insert(first_ele);
        // }
        uni_iters.emplace_back(++(array.cbegin()));
    }

    // start union
    auto in_iter = intersect.cbegin();
    while (!topns.empty()) {
        auto top_ele = topns.begin();
        auto& union_ele = top_ele->first;
        while (in_iter != intersect.cend() && *in_iter < union_ele) in_iter++;
        if (in_iter == intersect.cend()) break;
        // *in_iter&union_ele compare (len(union(unions))+len(intersect)) times
        if (*in_iter == union_ele) {
            in_iter++;
            for (auto& new_index : top_ele->second) {
                if (uni_iters[new_index] != unions[new_index].cend()) {
                    auto new_ele = *(uni_iters[new_index]++);
                    topns[new_ele].emplace_back(new_index);
                    // auto topn = topns.find(new_ele);
                    // if (topn != topns.end()) {
                    //     topn->second.emplace_back(new_index);
                    // } else {
                    //     pair<int, vector<int>> new_top({new_ele, {new_index}});
                    //     new_top.second.reserve(unions.size());
                    //     topns.insert(new_top);
                    // }
                }
                results[new_index].emplace(union_ele);
            }
        }
        else {
            for (auto& new_index : top_ele->second) {
                if (uni_iters[new_index] != unions[new_index].cend()) {
                    auto new_ele = *(uni_iters[new_index]++);
                    topns[new_ele].emplace_back(new_index);
                    // auto topn = topns.find(new_ele);
                    // if (topn != topns.end()) {
                    //     topn->second.emplace_back(new_index);
                    // } else {
                    //     pair<int, vector<int>> new_top({new_ele, {new_index}});
                    //     new_top.second.reserve(unions.size());
                    //     topns.insert(new_top);
                    // }
                }
            }
        }
        topns.erase(top_ele);
    }

    return results;
}

vector<set<int>>
union_concur_intersect_queue(const vector<set<int>>& unions, const set<int>& intersect) {
    vector<set<int>> results;
    results.resize(unions.size(), {});
    priority_queue<pair<int, int>, vector<pair<int, int>>, Compare> topns;
    pair<int, bool> last_value = { -1, false };

    vector<set<int>::iterator> uni_iters;
    for (int i = 0; i < unions.size(); i++) {
        auto& array = unions[i];
        int first_ele = *(array.begin());
        topns.push(make_pair(first_ele, i));
        uni_iters.emplace_back(++(array.cbegin()));
    }

    auto in_iter = intersect.cbegin();
    while (!topns.empty()) {
        auto top_ele = topns.top();
        if (top_ele.first == last_value.first) {
            auto& new_index = top_ele.second;
            if (last_value.second == true)
                results[new_index].emplace(top_ele.first);
            if (uni_iters[new_index] != unions[new_index].cend()) {
                auto new_ele = *(uni_iters[new_index]++);
                topns.push(make_pair(new_ele, new_index));
            }
            topns.pop();
            continue;
        }
        auto& union_ele = top_ele.first;
        last_value.first = union_ele;
        while (in_iter != intersect.cend() && *in_iter < union_ele) in_iter++;
        if (in_iter == intersect.cend()) break;
        auto& new_index = top_ele.second;
        if (uni_iters[new_index] != unions[new_index].cend()) {
            auto new_ele = *(uni_iters[new_index]++);
            topns.push(make_pair(new_ele, new_index));
        }
        if (*in_iter == union_ele) {
            in_iter++;
            results[new_index].emplace(union_ele);
            last_value.second = true;
        } else {
            last_value.second = false;
        }
        topns.pop();
    }

    return results;
}

void
union_concur_intersect_queue_vector(const vector<set<int>>& unions, const set<int>& intersect, vector<vector<int>>&results) {
    results.resize(unions.size(), {});
    priority_queue<pair<int, int>, vector<pair<int, int>>, Compare> topns;
    pair<int, bool> last_value = { -1, false };

    vector<set<int>::iterator> uni_iters;
    for (int i = 0; i < unions.size(); i++) {
        auto& array = unions[i];
        int first_ele = *(array.begin());
        topns.push(make_pair(first_ele, i));
        uni_iters.emplace_back(++(array.cbegin()));
    }

    auto in_iter = intersect.cbegin();
    while (!topns.empty()) {
        auto top_ele = topns.top();
        topns.pop();
        if (top_ele.first == last_value.first) {
            auto& new_index = top_ele.second;
            if (last_value.second == true)
                results[new_index].emplace_back(top_ele.first);
            if (uni_iters[new_index] != unions[new_index].cend()) {
                auto new_ele = *(uni_iters[new_index]++);
                topns.push(make_pair(new_ele, new_index));
            }
            continue;
        }
        auto& union_ele = top_ele.first;
        last_value.first = union_ele;
        while (in_iter != intersect.cend() && *in_iter < union_ele) in_iter++;
        if (in_iter == intersect.cend()) break;
        auto& new_index = top_ele.second;
        if (uni_iters[new_index] != unions[new_index].cend()) {
            auto new_ele = *(uni_iters[new_index]++);
            topns.push(make_pair(new_ele, new_index));
        }
        if (*in_iter == union_ele) {
            in_iter++;
            results[new_index].emplace_back(union_ele);
            last_value.second = true;
        } else {
            last_value.second = false;
        }
    }
}

#endif
