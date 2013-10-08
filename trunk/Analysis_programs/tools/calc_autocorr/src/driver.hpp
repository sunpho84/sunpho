#ifndef _DRIVER_H
#define _DRIVER_H

class driver_t
{
public:
  void *scanner;
  FILE *fin;
  driver_t(const char *path);
  virtual ~driver_t(){destroy_scanner();}
protected:
  void init_scanner();
  void destroy_scanner();
};

int parser_parse(driver_t *driver);
void read_from_file(const char *path);

#endif
