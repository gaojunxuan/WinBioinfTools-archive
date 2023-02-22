
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/uio.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>

#include "memmap.hpp"

using namespace std;

/*EE
  The following function returns a memory map for a given filename, or
  \texttt{NULL} if something went wrong.
*/

memmap::memmap(){
}

memmap::~memmap(){
}


void* memmap::creatememorymap(char *file,int line,char *filename, BOOL writemap,Uint *numofbytes)
{
  int fd;

//  DEBUG3(2,"\n# creatememorymap(file=%s,line=%d) for file %s:\n",file,line,filename);
  fd = simplefileOpen(filename,numofbytes);
  if(fd < 0)
  {
    return NULL;
  }
  return creatememorymapforfiledesc(file,line,fd,writemap,*numofbytes);
}


int memmap::deletememorymap(char *file,int line,void *mappedfile)
{
  int fd;
  fd=file_handle;

  if(mappedfile == NULL)
  {
    printf("%s: l. %d: deletememorymap: mappedfile is NULL",file,line);
    return -1;
  }
/*
  for(fd=0; fd<MAXMAPPEDFILES; fd++)
  {
    if(memoryptr[fd] == mappedfile)
    {
      break;
    }
  }
  if(fd == MAXMAPPEDFILES)
  {
    ERROR2("%s: l. %d: deletememorymap: cannot find filedescriptor for given address",
           file,line);
    return -2;
  }
*/
  if(munmap((caddr_t) memoryptr,(size_t) mappedbytes) != 0)
  {
    printf("Delete memory map failed \n");
    return -3;
  }
/*
  DEBUG4(2,"# file \"%s\", line %d: deletememorymap:fd=%d (%u) bytes:\n",
          file,line,fd,mappedbytes[fd]);
  DEBUG2(2,"# mapped in file \"%s\", line %d\n",filemapped[fd],
                                                linemapped[fd]);
*/
  memoryptr = NULL;
//  mmsubtractspace(mappedbytes[fd]);
  mappedbytes = 0;
  close(fd);
  return 0;
}

/*EE
  The following function creates a memory map for a given file
  descriptor \texttt{fd}. \texttt{writemap} is true iff the map
  should be writable, and \texttt{numofbytes} is the 
  size of the file to be mapped in bytes.
*/

void* memmap::creatememorymapforfiledesc(char *file,int line,int fd,BOOL writemap,Uint numofbytes)
{
    
//  DEBUG1(2,"# creatememorymapforfiledesc %d\n",fd);
/* 
 if(fd < 0)
  {
    printf("creatememorymap: filedescriptor %d negative",fd);
    return NULL;
  }
  if(fd >= MAXMAPPEDFILES)
  {
    printf("creatememorymap: filedescriptor %d is too large",fd);
    return NULL;
  }
  if(memoryptr[fd] != NULL)
  {
    printf("creatememorymap: filedescriptor %d already in use",fd);
    return NULL;
  }

  mmaddspace(numofbytes);

  //DEBUG2(2,"# memorymap:fd=%d: %u bytes\n",fd,numofbytes);
  mappedbytes[fd] = numofbytes;
*/
    mappedbytes = numofbytes;
    memoryptr= (void *) mmap(0,(size_t) numofbytes,
                                writemap ? (PROT_READ | PROT_WRITE) : PROT_READ,
                                MAP_PRIVATE,
                                fd,(off_t) 0);
  if(memoryptr == (void *) MAP_FAILED)
  {
    printf("memorymapping for filedescriptor %d failed",fd);
    return NULL;
  }
//  filemapped[fd] = file;
//  linemapped[fd] = line;
  file_handle=fd;
  return memoryptr;
}


int memmap:: simplefileOpen(char *filename,Uint *numofbytes)
{
  int fd;
  struct stat buf;

  if((fd = open(filename,O_RDONLY)) == -1)
  {
     printf("cannot open \"%s\"",filename);
     return -1;
  }
  if(fstat(fd,&buf) == -1)
  {
     printf("cannot access file status for \"%s\"",filename);
     return -2;
  }
  *numofbytes = (Uint) buf.st_size;
  return fd;
}
