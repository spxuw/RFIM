#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
using namespace std;


void MakeDirectories()
{
    struct stat s;
    if (-1==stat("data",&s))
    {
	int  err = system("mkdir data");
	if (err == -1) 
	    cout << "Failed to create directory data " <<endl;
    }
}


void MakeDirectories(double m)
{
    struct stat s;
 
    char temp [256];
    sprintf(temp, "data/m%.3lf/", m);
    char cmd [256];
    sprintf(cmd, "mkdir %s", temp);

   if (-1==stat(temp,&s))
    {
	int  err = system(cmd);
	if (err == -1) 
	    cout << "Failed to create dir " << temp << endl;
    }
}




// for real data analysis
void MakeDirectory(string path, int type)
{
  struct stat s;
  
  char temp [256];
  sprintf(temp, "%s/rand%d", path.c_str(), type);
  char cmd [256];
  sprintf(cmd, "mkdir %s", temp);
  if (-1==stat(temp,&s)) {
    int  err = system(cmd);
    if (err == -1) 
      cout << "Failed to create dir " << temp << endl;
  }

}


// for real data analysis
void MakeDirectory(string path, double p)
{
  struct stat s;
  
  char temp [256];
  sprintf(temp, "%s/p%.3lf", path.c_str(), p);
  char cmd [256];
  sprintf(cmd, "mkdir %s", temp);
  if (-1==stat(temp,&s)) {
    int  err = system(cmd);
    if (err == -1) 
      cout << "Failed to create dir " << temp << endl;
  }

}


// for real data analysis
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


// for real data analysis
void MakeDirectories(string file)
{
  struct stat s;
  /*
  if (-1==stat("data",&s))  {
    int  err = system("mkdir data");
    if (err == -1) 
      cout << "Failed to create directory data " <<endl;
  }
  */
  
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



// for synthetic data analysis
void MakeDirectories(char* file)
{
    struct stat s;

    
    if (-1==stat("data",&s))  {
      int  err = system("mkdir data");
      if (err == -1) 
	cout << "Failed to create directory data " <<endl;
    }
    
    char temp [256];
    sprintf(temp, "data/%s", file);
    char cmd [256];
    sprintf(cmd, "mkdir %s", temp);

   if (-1==stat(temp,&s))
    {
	int  err = system(cmd);
	if (err == -1) 
	    cout << "Failed to create dir " << temp << endl;
    }
}


// for synthetic data analysis
void MakeDirectories(char* file, int seed)
{
    struct stat s;

    char temp [256];
    sprintf(temp, "data/%s/seed%d/", file,seed);
    char cmd [256];
    sprintf(cmd, "mkdir %s", temp);

   if (-1==stat(temp,&s))
    {
	int  err = system(cmd);
	if (err == -1) 
	    cout << "Failed to create dir " << temp << endl;
    }
}


