#ifdef HAVE_CONFIG_H
 #include "config.hpp"
#endif

//check the endianness of the machine
bool get_little_endianness()
{
  int little_endian=1;
  little_endian=(int)(*(char*)(&little_endian));
  return (little_endian==1);
}
