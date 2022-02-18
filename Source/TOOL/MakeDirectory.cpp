#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
using namespace std;






void MakeDirectory()
{
    struct stat s;
    if (-1==stat("data",&s))
    {
	int  err = system("mkdir data");
	if (err == -1) 
	    cout << "Failed to create directory data " <<endl;
    }
    if (-1==stat("data/GSlayer",&s))
    {
	int  err = system("mkdir data/GSlayer");
	if (err == -1) 
	    cout << "Failed to create dir data/GSlayer" << endl;
    }

}






void MakeDirectory_Default()
{
    struct stat s;
    if (-1==stat("data",&s))
    {
	int  err = system("mkdir data");
	if (err == -1) 
	    cout << "Failed to create directory data " <<endl;
    }
}





void MakeDirectory(int D, int L, double R, int seed)
{
    struct stat s;
    if (-1==stat("data",&s))
    {
	int  err = system("mkdir data");
	if (err == -1) 
	    cout << "Failed to create directory data " <<endl;
    }

    char jobtmp[256]; 
    sprintf(jobtmp,"./data/GS-D%d-L%d-R%.3lf-seed%d",D,L,R,seed);
    char cmd[256];
    sprintf(cmd,"mkdir %s",jobtmp);
    if (-1==stat(jobtmp,&s))
    {
	int  err = system(cmd);
	if (err == -1) 
	    cout << "Failed to create dir" << jobtmp << endl;
    }
}




void MakeDirectory(int N, double c, double R, int seed)
{
    struct stat s;
    if (-1==stat("data",&s))
    {
	int  err = system("mkdir data");
	if (err == -1) 
	    cout << "Failed to create directory data " <<endl;
    }

    char jobtmp[256]; 
    sprintf(jobtmp,"./data/GS-N%d-c%e-R%.3lf-seed%d",N,c,R,seed);
    char cmd[256];
    sprintf(cmd,"mkdir %s",jobtmp);
    if (-1==stat(jobtmp,&s))
    {
	int  err = system(cmd);
	if (err == -1) 
	    cout << "Failed to create dir" << jobtmp << endl;
    }
}







void MakeDirectory(string file)
{
  struct stat s;

  char temp [256];
  sprintf(temp, "data-%s/", file.c_str());
  char cmd [256];
  sprintf(cmd, "mkdir %s", temp);
  if (-1==stat(temp,&s)) {
    int  err = system(cmd);
    if (err == -1) 
      cout << "Failed to create dir " << temp << endl;
  }

}



void MakeDirectories(string file)
{
  struct stat s;
 
  
  char temp [256];
  sprintf(temp, "data-%s/", file.c_str());
  char cmd [256];
  sprintf(cmd, "mkdir %s", temp);
  if (-1==stat(temp,&s)) {
    int  err = system(cmd);
    if (err == -1) 
      cout << "Failed to create dir " << temp << endl;
  }


  sprintf(temp, "data-%s/Real/", file.c_str());
  sprintf(cmd, "mkdir %s", temp);
  if (-1==stat(temp,&s)) {
    int  err = system(cmd);
    if (err == -1) 
      cout << "Failed to create dir " << temp << endl;
  }
  

  sprintf(temp, "data-%s/Rand0/", file.c_str());
  sprintf(cmd, "mkdir %s", temp);
  if (-1==stat(temp,&s)) {
    int  err = system(cmd);
    if (err == -1) 
      cout << "Failed to create dir " << temp << endl;
  }


  sprintf(temp, "data-%s/Rand1/", file.c_str());
  sprintf(cmd, "mkdir %s", temp);
  if (-1==stat(temp,&s)) {
    int  err = system(cmd);
    if (err == -1) 
      cout << "Failed to create dir " << temp << endl;
  }


  sprintf(temp, "data-%s/Rand2/", file.c_str());
  sprintf(cmd, "mkdir %s", temp);
  if (-1==stat(temp,&s)) {
    int  err = system(cmd);
    if (err == -1) 
      cout << "Failed to create dir " << temp << endl;
  }


}

