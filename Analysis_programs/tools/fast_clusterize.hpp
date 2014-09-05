#include <array>
#include <cstdlib>
#include <iostream>
#include <fstream>

template<class T> T revert_endianness(T a)
{
  union rev_t
  {
    T num;
    char car[sizeof(T)];
  };
  
  rev_t in;
  rev_t out;
  
  in.num=a;
  for(int i=0;i<sizeof(T);i++) out.car[i]=in.car[(sizeof(T)-i)%sizeof(T)];
  
  return out.num;
}

void crash(const char *temp,...)
{
  char buffer[1024];
  va_list args;

  va_start(args,temp);
  vsprintf(buffer,temp,args);
  va_end(args);

  std::cerr<<"ERROR: "<<buffer<<std::endl;
  exit(1);
}

template <int njacks> class jack_t : public std::array<double,njacks+1>
{
 private:
 public:
  jack_t();
  void clusterize(int nconfs);
};

template<int njacks> jack_t<njacks>::jack_t()
{for(int ijack=0;ijack<1;ijack++) (*this).operator[](ijack)=0;}

template<int njacks> void jack_t<njacks>::clusterize(int nconfs)
{
  if(nconfs%njacks) crash("%d/%d!=0",nconfs,njacks);
  int clust_size=nconfs/njacks;

  (*this)[njacks]=0;
  for(int ijack=0;ijack<njacks;ijack++) (*this)[njacks]+=(*this)[ijack];
  for(int ijack=0;ijack<njacks;ijack++) (*this)[ijack]=((*this)[njacks]-(*this)[ijack])/(nconfs-clust_size);
  (*this)[njacks]/=nconfs;
}

class sm_ifstream : public std::ifstream
{
private:
  bool end_flag;
  sm_ifstream();
public:
  void reverse_endianness(){end_flag=!end_flag;}
  sm_ifstream(std::string path_in)
  {
    end_flag=0;
    open(path_in);
    if(!good()) crash("opening %s",path_in.c_str());
  }
  template<class T> T read()
  {
    T out;
    
    if(!std::ifstream::read((char*)&out,sizeof(T))) crash("reading type of size %u",sizeof(T));

    if(end_flag) return revert_endianness(out);
    else return out;
  }
  streampos get_size()
  {
    streampos ori=tellg();
    seekg(0,ios_base::end);
    streampos end=tellg();
    seekg(ori);
    
    return end-ori;
  }
};

class sm_ofstream : public std::ofstream
{
private:
  int end_flag;
  sm_ofstream();
public:
  sm_ofstream(std::string path_out)
  {
    end_flag=0;
    open(path_out);
    if(!good()) crash("opening %s",path_out.c_str());
  }
  template<class T> void write(T &in)
  {
    T *p;
    if(end_flag)
      {
	p=new T;
	(*p)=revert_endianness(in);
      }
    else p=&in;
    
    if(!std::ofstream::write((char*)p,sizeof(T))) crash("writing type of size %u",sizeof(T));
    if(end_flag) delete p;
  }

};
