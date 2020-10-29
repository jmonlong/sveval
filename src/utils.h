#ifndef UTILS_H
#define UTILS_H 1

// split a string
std::vector<std::string> split_str(std::string in_str, std::string delim=std::string("\t"));

// reverse complement
std::string rev_comp(std::string seq);

#endif // UTILS_H

