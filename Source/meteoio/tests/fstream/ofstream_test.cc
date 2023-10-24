#include <stdlib.h>
#include <meteoio/MeteoIO.h>
#include <sys/stat.h>
#include <fstream>
#include <meteoio/FStream.h>


using namespace mio; //The MeteoIO namespace is called mio

int main() {
    ofilestream tmp;
    bool dir_status, filestatus=true;
    if (tmp.getDefault()) {
        dir_status = false;
        filestatus = false;
        std::cerr<< "Testing the generation of directories from path"<<std::endl;
        std::string filename1 = "./this/is/a/test.txt";
        std::string filename2 = "./this/is/a/second_test.txt";
        ofilestream ofs1(filename1.c_str());
        ofilestream ofs2(filename2.c_str());
        std::string instr ="This is a test";
        ofs1 << instr;
        ofs2 << instr;
        ofs1.close();
        ofs2.close();
        struct stat sb;

        if (stat(mio::FileUtils::getPath(filename1,false).c_str(), &sb) == 0) {
            dir_status = true;
        } else {
            throw IOException("Output directory was not created.");
        }

        std::ifstream infile1(filename1);
        std::ifstream infile2(filename2);
        std::string item;
        while (getline(infile1, item)) {
            if (item != instr) {
                throw IOException("File was not created, or written to!");
            } else {
                filestatus = true;
            }
        }
        while (getline(infile2, item)) {
            if (item != instr) {
                throw IOException("File was not created, or written to!");
            } else {
                filestatus = true;
            }
        }
        std::string win_path = "C:\\Users\\Thomas\\Documents\\MeteoIO\\tests\\fstream\\test.txt";
        std::string win_path2 = "C:/Users/Thomas/Documents/MeteoIO/tests/fstream/test.txt";
        try
        {
            ofilestream win(win_path);
            std::cerr << "Creating Windows path type directories should not be possible";
            exit(1);
        }
        catch(const std::exception& e)
        {
        }
        try
        {
            ofilestream win2(win_path2);
            std::cerr << "Creating Windows path type directories should not be possible";
            exit(1);
        }
        catch(const std::exception& e)
        {
        }

        std::string abs_path3 = "/home/leibersp/Documents/test.txt"; //should show a warning but work
        ofilestream abs_fail(abs_path3);
        abs_fail << "hiii";
        abs_fail.close();
    
        system("rm -rf this");
        system("rm -rf Documents");

        std::cerr << "Everything works as expected"<< std::endl;
        std::cerr << "----------------------------" << std::endl;

#ifdef LIMIT_WRITE_ACCESS
        std::cerr << "Testing to write a not allowed directory with limited writing access"<<std::endl;
        std::string filename3 ="../../this_shouldnt_be_outside_of_fstream/works.txt";
        ofilestream file3(filename3);
        file3 << "hiii";
        file3.close();
        if (FileUtils::directoryExists("../../this_shouldnt_be_outside_of_fstream")) {
            std::cerr << "limiting write acces did not work"<<std::endl;
        } else if (FileUtils::directoryExists("this_shouldnt_be_outside_of_fstream")) {
            std::cerr << "limiting write acess works perfectly"<< std::endl;
        } else {
            std::cerr << "directory was not created at all in write acces limitation scope"<< std::endl;
        }
        system ("rm -rf this_shouldnt_be_outside_of_fstream");
        std::cerr << "----------------------------" << std::endl;
#endif
    } else {
        std::cerr << "Writing of directories not enabled"<<std::endl;
        std::cerr << "Testing the generation of files"<< std::endl;
        std::string filename = "trial.txt";
        ofilestream file(filename);
        file << "nice" << std::endl;
        file.close();
        std::ifstream infile(filename);
        std::string item;
        while (getline(infile,item)) {
            if (item != "nice"){
                std::cerr << "File Was not created or written to";
                exit(1);
            }
        }
        system("rm trial.txt");
    }

    std::cerr << "checking the reading of ini file"<<std::endl;
    mio::Config cfg("io.ini");

    if (cfg.get("WRITE_DIRECTORIES","Output",true)!=false) throw IOException("Ini file was not read properly");
    IOManager mng(cfg);
    ofilestream tmp2("hi");
    system("rm hi");
    if (tmp2.getDefault() != false) throw IOException("Default was not set");
#ifdef LIMIT_WRITE_ACCESS
    std::string abs_path1 = "/home/leibersp/Documents/test.txt";// should not work
    try
    {
        ofilestream fail(abs_path1);
        std::cerr<< "Default |" << fail.getDefault() << std::endl;
        if (fail.getDefault()) std::cerr<<"Default was not kept"<<std::endl;
        std::cerr << "accessing absolute invalid path, without writing directories should not be possible"<<std::endl;
        exit(1);
    }
    catch(const std::exception& e)
    {
        std::cerr << "accessing absolute invalid path, without writing directories works as expected"<< std::endl;
    }
    
    std::string abs_path2 = "/home/leibersp/slf/meteoio/tests/fstream/test.txt";
    ofilestream works(abs_path2);
    if (works.getDefault()) throw IOException("Default was not kept");
    works << "hiii";
    works.close();
    system("rm test.txt");
#endif

    std::cerr << "Checking for calls to std::ofstream outside of wrapper"<< std::endl;
    system("rgrep --exclude=\"*.cc.o\" \"ofstream\" ../../meteoio | grep -vE \"^../../meteoio/FStream.\" | grep -vE \"ofstream::app\" |grep -vE \"ofstream::out\">tmp_ofstream_instances.txt");
    std::ifstream ofstream_instances("tmp_ofstream_instances.txt");
    std::string item;
    while (getline(ofstream_instances, item,':')) {
        if (item.empty()) continue;
        else {
            std::cerr <<"Instance of std::ofstream found in file " <<item;
            exit(1);
        }
    }
    system("rm tmp_ofstream_instances.txt");
    std::cerr << "Nothing found"<< std::endl;

    return !filestatus || !dir_status;
}
