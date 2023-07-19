#include <iostream>
#include <fstream>
#include <string>


int main(int argc, char *argv[])
{
  if (argc != 2) { std::cout << "Need one file name as argument!" << std::endl; return 0; }

  std::string inputName  = argv[1];
  std::string outputName = "NEW_"+inputName;

  std::ifstream inputFile(inputName.c_str());
  std::ofstream outputFile(outputName.c_str());

  std::string tempLine;
  std::string tempPath;

  while (inputFile.good() && inputFile >> tempLine)
    {
      if (isdigit(tempLine[0])) continue;
      else outputFile << tempLine << std::endl;
    }

  inputFile.close();
  outputFile.close();
  return 0;
}
