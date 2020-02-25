import sys
import configparser
import numpy as np

class IniFile:
    
    def __init__(self,fname):
        self.ini = configparser.ConfigParser()
        self.ini.read(fname)
        
    def Dump(self):
        print("# inifile:")
        sections = self.ini.sections()
        for section in sections:
            print(' ['+section+']')
            for key in self.ini[section]:
                print(" ",key,"=",self.ini[section][key])
        print("")
          
    def ReadBoolean(self,section,key):
        return self.ini[section].getboolean(key)
    
    def ReadFloat(self,section,key):
        return float(self.ini[section][key])
    
    def ReadFloatArray(self,section,key):
        s = self.ini[section][key]
        return np.fromstring(s,dtype=float,sep=",")

    def ReadInt(self,section,key):
        return int(self.ini[section][key])
    
    def ReadString(self,section,key):
        return self.ini[section][key]
    
def main():
    args = sys.argv
    if(len(args)<2):
        print("error: the number of input parameters is not correct; input command must be")
        print(">$ python inifiles.py input_file_name")
        sys.exit(1)
        
    Ini = IniFile(args[1])
    Ini.Dump("test.out")

if __name__=='__main__':
    main()
