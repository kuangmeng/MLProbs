#include "qscore.h"

int main(int argc, char *argv[])
	{
	if (argc < 2)
		{
		Usage();
		exit(0);
		}

	ParseOptions(argc - 1, argv + 1);
	if (FlagOpt("version"))
		{
		fprintf(stderr, "qscore v2.1\n");
		exit(0);
		}

	if (ValueOpt("sab_test"))
		{
		SAB();
		return 0;
		}

	QScore();
	return 0;
	}
