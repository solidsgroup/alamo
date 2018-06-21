#include "FileNameParse.H"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>

std::time_t t = 0;

int ReplaceAll(std::string &str, const std::string before, const std::string after)
{
  while(str.find(before) != std::string::npos)
    str.replace(str.find(before), before.size(), after);
  return 0;
}

void IO::FileNameParse(std::string &filename)
{
	if (!t) t = std::time(0);
	std::tm * now = std::localtime(&t);
	int year = now->tm_year+1900;
	int month = now->tm_mon+1;
	int day = now->tm_mday;
	int hour = now->tm_hour;
	int minute = now->tm_min;
	int second = now->tm_sec;
	
	std::stringstream ss;

	// _Y: year
	ss.str("");
	ss << year;
	std::string _Y = ss.str();
	
	// _m: month (01..12)
	ss.str("");
	ss << std::setfill('0') << std::setw(2) << month;
	std::string _m = ss.str();
	
	// _d: day of month (01..31)
	ss.str("");
	ss << std::setfill('0') << std::setw(2) << day;
	std::string _d = ss.str();
	
	// _H: hour (00..23)
	ss.str("");
	ss << std::setfill('0') << std::setw(2) << hour;
	std::string _H = ss.str();
	
	// _M: minute (00..59)
	ss.str("");
	ss << std::setfill('0') << std::setw(2) << minute;
	std::string _M = ss.str();

	// _S: second (00..59)
	ss.str("");
	ss << std::setfill('0') << std::setw(2) << second;
	std::string _S = ss.str();

	ReplaceAll(filename,"%Y",_Y);
	ReplaceAll(filename,"%m",_m);
	ReplaceAll(filename,"%d",_d);
	ReplaceAll(filename,"%H",_H);
	ReplaceAll(filename,"%M",_M);
	ReplaceAll(filename,"%S",_S);
}
