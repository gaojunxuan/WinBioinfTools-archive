#ifndef MEMMAP_H
#define MEMMAP_H

//#ifndef BOOL char
#define BOOL char
//#endif

//#ifndef True
#define True 1
//#endif

//#ifndef False
#define False 0
//#endif

//#ifndef Uint
#define Uint long
//#endif

class memmap{

public:
    memmap();

    ~memmap();
public:
   void *memoryptr;
   long mappedbytes;
    int file_handle;
    void* creatememorymap(char *file,int line,char *filename, BOOL writemap,Uint *numofbytes);
    void* creatememorymapforfiledesc(char *file,int line,int fd, BOOL writemap,Uint numofbytes);
    int deletememorymap(char *file,int line,void *mappedfile);
    int simplefileOpen(char *filename,Uint *numofbytes);
};

#define CONSTRUCTMEMMAP()\
        memmap memmap_obj

#define CREATEMEMORYMAP(F,WM,NB)\
        memmap_obj.creatememorymap(__FILE__,__LINE__,F,WM,NB)

#define CREATEMEMORYMAPFORFILEDESC(FD,WM,NB)\
        creatememorymapforfiledesc(__FILE__,__LINE__,FD,WM,NB)

#define DELETEMEMORYMAP(MF)\
        memmap_obj.deletememorymap(__FILE__,__LINE__,(void *) MF)

//\Ignore{

#endif

//}
