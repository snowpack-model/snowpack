#include <stdlib.h>
#include <meteoio/MeteoIO.h>
#include <sys/stat.h>
#include <fstream>
#include <meteoio/FileUtils.h>
#include <meteoio/FStream.h>


using namespace mio; //The MeteoIO namespace is called mio

static void testRelativePath(bool& dir_status, bool& filestatus)
{
	dir_status = false;
	filestatus = false;

	std::cerr << "Testing the generation of directories from path" << std::endl;
	static const std::string filename1 = "./this/is/a/test.txt";
	static const std::string filename2 = "./this/is/a/second_test.txt";
	ofilestream ofs1(filename1.c_str());
	ofilestream ofs2(filename2.c_str());
	std::string instr = "This is a test";
	ofs1 << instr;
	ofs2 << instr;
	ofs1.close();
	ofs2.close();

	struct stat sb;
	if (stat(mio::FileUtils::getPath(filename1, false).c_str(), &sb) == 0) {
		dir_status = true;
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

	infile1.close();
	infile2.close();
}

static void testWinPath()
{
	static const std::string win_path = "C:\\Users\\Thomas\\Documents\\MeteoIO\\tests\\fstream\\test.txt";
	static const std::string win_path2 = "C:/Users/Thomas/Documents/MeteoIO/tests/fstream/test.txt";
	
	try {
		ofilestream win(win_path);
		std::cerr << "Creating Windows path type directories should not be possible";
		exit(1);
	} catch(const std::exception&) { }
	
	try {
		ofilestream win2(win_path2);
		std::cerr << "Creating Windows path type directories should not be possible";
		exit(1);
	} catch(const std::exception&) { }
}

static void testAbsPath()
{
	static const std::string abs_path3 = mio::FileUtils::cleanPath("./", true) + "rand3/test.txt";
	mio::FileUtils::createDirectories(abs_path3);

	try {
		ofilestream abs_fail(abs_path3.c_str());
		abs_fail << "hiii";
		abs_fail.close();
	}
	catch (const std::exception&) {
		// Handle the exception if needed
	}

	system(("rm -rf " + abs_path3).c_str());
}

#ifdef LIMIT_WRITE_ACCESS
static void testLimitedWriteAccess()
{
	std::cerr << "Testing to write a not allowed directory with limited writing access" << std::endl;
	static const std::string dirname3 = "this_shouldnt_be_outside_of_fstream";
	static const std::string filename3 = "../../" + dirname3 + "/works.txt";
	ofilestream file3(filename3);
	file3 << "hiii";
	file3.close();
	if (FileUtils::directoryExists("../../" + dirname3))
	{
		std::cerr << "Limiting write access did not work" << std::endl;
	}
	else if (FileUtils::directoryExists(dirname3))
	{
		std::cerr << "Limiting write access works perfectly" << std::endl;
	}
	else
	{
		std::cerr << "Directory was not created at all in write access limitation scope" << std::endl;
	}
	system(("rm -rf " + dirname3).c_str());
}

static void testLimitedAccessWithoutCreatingDirectories() {
	const std::string abs_path11 = mio::FileUtils::cleanPath("./", true) + "rand/test.txt"; // should not work
	try {
		ofilestream fail(abs_path1);
		std::cerr << "Default |" << fail.getDefault() << std::endl;
		if (fail.getDefault())
			std::cerr << "Default was not kept" << std::endl;
		std::cerr << "Accessing absolute invalid path, without writing directories should not be possible" << std::endl;
		fail.close();
		exit(1);
	} catch (const std::exception &e) {
		std::cerr << e.print() << "\n";
		std::cerr << "Accessing absolute invalid path, without writing directories works as expected" << std::endl;
	}
	system("rm -rf " + abs_path11);

	const std::string abs_path2 = mio::FileUtils::cleanPath("./", true) + "rand2/test.txt";
	mio::FileUtils::createDirectories(abs_path2);
	ofilestream works(abs_path2);
	if (works.getDefault())
		throw IOException("Default was not kept");
	works << "hiii";
	works.close();
	system("rm -rf " + abs_path2);
	std:: cerr << "Accessing absolute valid path, without writing directories works as expected" << std::endl;
}
#endif

static void testFileWrite()
{
	std::cerr << "Writing of directories not enabled" << std::endl;
	std::cerr << "Testing the generation of files" << std::endl;
	static const std::string filename = "trial.txt";
	ofilestream file(filename);
	file << "nice" << std::endl;
	file.close();
	std::ifstream infile(filename);
	std::string item;
	while (getline(infile, item))
	{
		if (item != "nice")
		{
			std::cerr << "File was not created or written to" << std::endl;
			exit(1);
		}
	}
	infile.close();
	std::remove(filename.c_str());
}
static void searchForRawCalls()
{
	std::cerr << "Checking for calls to std::ofstream outside of wrapper" << std::endl;
	// Exclude occurrences of "ofstream" that are commented out
	system("rgrep --exclude=\"*.cc.o\" \"^[^\\*\\/]*(ofstream)\" ../../meteoio | grep -vE \"^../../meteoio/FStream.\" | grep -vE \"ofstream::app\" | grep -vE \"ofstream::out\">tmp_ofstream_instances.txt");
	std::ifstream ofstream_instances("tmp_ofstream_instances.txt");
	std::string item;
	while (getline(ofstream_instances, item, ':')) {
		if (item.empty())
			continue;
		else {
			std::cerr << "Instance of std::ofstream found in file " << item;
			ofstream_instances.close();
			exit(1);
		}
	}
	ofstream_instances.close();
	std::remove("tmp_ofstream_instances.txt");
	std::cerr << "Nothing found" << std::endl;
}

int main() {
	bool dir_status = true, filestatus = true;
	if (ofilestream().getDefault()) {
		testRelativePath(dir_status, filestatus);
		testWinPath();
		testAbsPath();

		system("rm -rf ./this");
		system("rm -rf ./Documents");

		std::cerr << "Everything works as expected" << std::endl;
		std::cerr << "----------------------------" << std::endl;

#ifdef LIMIT_WRITE_ACCESS
		testLimitedWriteAccess();
		std::cerr << "----------------------------" << std::endl;
#endif
	} else {
		testFileWrite();
	}

	std::cerr << "Checking the reading of ini file" << std::endl;
	mio::Config cfg("io.ini");

	if (cfg.get("WRITE_DIRECTORIES", "Output", true) != false)
		throw IOException("Ini file was not read properly");

	IOManager mng(cfg);
	ofilestream tmp2("hi");
	system("rm hi");
	if (tmp2.getDefault() != false)
		throw IOException("Default was not set");

#ifdef LIMIT_WRITE_ACCESS
	testLimitedAccessWithoutCreatingDirectories();
#endif

	searchForRawCalls();

	return !filestatus || !dir_status;
}
