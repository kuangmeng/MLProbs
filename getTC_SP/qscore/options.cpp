#include "qscore.h"

struct VALUE_OPT
	{
	const char *m_pstrName;
	const char *m_pstrValue;
	};

struct FLAG_OPT
	{
	const char *m_pstrName;
	bool m_bSet;
	};

static VALUE_OPT ValueOpts[] =
	{
	{ "test",					0, },
	{ "ref",					0, },
	{ "sab_test",				0, },
	{ "sab_ref",				0, },
	};
static int ValueOptCount = sizeof(ValueOpts)/sizeof(ValueOpts[0]);

static FLAG_OPT FlagOpts[] =
	{
	{ "truncname",			false, },
	{ "ignoretestcase",		false, },
	{ "ignorerefcase",		false, },
	{ "quiet",				false, },
	{ "cline",				false, },
	{ "modeler",			false, },
	{ "slow",				false, },
	{ "version",			false, },
	{ "gapscore",			false, },
	{ "seqdiffwarn",		false, },
	{ "ignoremissingseqs",	false, },
	{ "perseq",				false, },
	{ "verbose",			false, },
	{ "stripx",				false, },
	{ "stripb",				false, },
	{ "stripz",				false, },
	};
static int FlagOptCount = sizeof(FlagOpts)/sizeof(FlagOpts[0]);

static bool TestSetFlagOpt(const char *Arg)
	{
	for (int i = 0; i < FlagOptCount; ++i)
		if (!strcmp(Arg, FlagOpts[i].m_pstrName))
			{
			FlagOpts[i].m_bSet = true;
			return true;
			}
	return false;
	}

static bool TestSetValueOpt(const char *Arg, const char *Value)
	{
	for (int i = 0; i < ValueOptCount; ++i)
		if (!strcmp(Arg, ValueOpts[i].m_pstrName))
			{
			if (0 == Value)
				{
				fprintf(stderr, "Option -%s must have value\n", Arg);
				exit(1);
				}
			ValueOpts[i].m_pstrValue = strdup(Value);
			return true;
			}
	return false;
	}

bool FlagOpt(const char *Name)
	{
	for (int i = 0; i < FlagOptCount; ++i)
		if (!strcmp(Name, FlagOpts[i].m_pstrName))
			return FlagOpts[i].m_bSet;
	Quit("FlagOpt(%s) invalid", Name);
	return false;
	}

const char *ValueOpt(const char *Name)
	{
	for (int i = 0; i < ValueOptCount; ++i)
		if (!strcmp(Name, ValueOpts[i].m_pstrName))
			return ValueOpts[i].m_pstrValue;
	Quit("ValueOpt(%s) invalid", Name);
	return 0;
	}

const char *RequiredValueOpt(const char *Name)
	{
	const char *s = ValueOpt(Name);
	if (0 == s)
		Quit("Required option -%s not specified", Name);
	return s;
	}

void ParseOptions(int argc, char *argv[])
	{
	for (int iArgIndex = 0; iArgIndex < argc; )
		{
		const char *Arg = argv[iArgIndex];
		if (Arg[0] != '-')
			Quit("Command-line option \"%s\" must start with '-'", Arg);
		const char *ArgName = Arg + 1;
		if (TestSetFlagOpt(ArgName))
			{
			++iArgIndex;
			continue;
			}
		
		char *Value = 0;
		if (iArgIndex < argc - 1)
			Value = argv[iArgIndex+1];
		if (TestSetValueOpt(ArgName, Value))
			{
			iArgIndex += 2;
			continue;
			}
		fprintf(stderr, "Invalid command line option \"%s\"\n", ArgName);
		Usage();
		exit(1);
		}
	}
