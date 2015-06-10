#include <sys/time.h>
#include <unistd.h>

struct myCPUtimer
{
	struct timeval m_start, m_end;
	float m_seconds, m_useconds;
};

float TMinit(struct myCPUtimer* timer)
{
	timer->m_end.tv_sec = 0;	
	timer->m_end.tv_usec = 0;	
	gettimeofday(&timer->m_start, NULL);
	timer->m_seconds = timer->m_start.tv_sec; 
	timer->m_useconds = timer->m_start.tv_usec; 
	return ((timer->m_seconds) * 1000 + timer->m_useconds/1000.0);
}
float TMstop(struct myCPUtimer* timer)
{
	if(timer->m_start.tv_sec == 0)
	{
		return -1;
	}
	gettimeofday(&timer->m_end, NULL);
	timer->m_seconds = timer->m_end.tv_sec - timer->m_start.tv_sec; 
	timer->m_useconds = timer->m_end.tv_usec - timer->m_start.tv_usec; 
	return ((timer->m_seconds) * 1000 + timer->m_useconds/1000.0);
}
float TMget_time(struct myCPUtimer* timer)
{
	if(timer->m_start.tv_sec == 0)
	{
		return -1;
	}
	if(timer->m_end.tv_sec == 0)
	{
		return TMstop(timer);
	}
	return ((timer->m_seconds) * 1000 + timer->m_useconds/1000.0);
}
float TMrestart(struct myCPUtimer* timer)
{
	if(timer->m_start.tv_sec == 0)
	{
		return -1;
	}
	gettimeofday(&timer->m_end, NULL);
	timer->m_seconds = timer->m_end.tv_sec - timer->m_start.tv_sec;
	timer->m_useconds = timer->m_end.tv_usec - timer->m_start.tv_usec;
	timer->m_start.tv_sec = timer->m_end.tv_sec;
	timer->m_start.tv_usec = timer->m_end.tv_usec;
	return ((timer->m_seconds) * 1000 + timer->m_useconds/1000.0);
}
