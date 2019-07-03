#include "FileNameParse.H"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>

#include "Util/Util.H"

std::time_t t = 0;

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
	Util::String::ReplaceAll(filename,"%Y",_Y);
	
	// _m: month (01..12)
	ss.str("");
	ss << std::setfill('0') << std::setw(2) << month;
	std::string _m = ss.str();
	Util::String::ReplaceAll(filename,"%m",_m);
	
	// _d: day of month (01..31)
	ss.str("");
	ss << std::setfill('0') << std::setw(2) << day;
	std::string _d = ss.str();
	Util::String::ReplaceAll(filename,"%d",_d);
	
	// _H: hour (00..23)
	ss.str("");
	ss << std::setfill('0') << std::setw(2) << hour;
	std::string _H = ss.str();
	Util::String::ReplaceAll(filename,"%H",_H);
	
	// _M: minute (00..59)
	ss.str("");
	ss << std::setfill('0') << std::setw(2) << minute;
	std::string _M = ss.str();
	Util::String::ReplaceAll(filename,"%M",_M);

	// _S: second (00..59)
	ss.str("");
	ss << std::setfill('0') << std::setw(2) << second;
	std::string _S = ss.str();
	Util::String::ReplaceAll(filename,"%S",_S);

	// _D: spatial dimension (1,2,3)
	ss.str("");
	ss << AMREX_SPACEDIM;
	std::string _D = ss.str();
	Util::String::ReplaceAll(filename,"%D",_D);
	
}
