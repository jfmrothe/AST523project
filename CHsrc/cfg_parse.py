#!/usr/bin/python
#cfg_parse.py use to parse all the configure parameters from .cfg file for 
#functions related to run lc function. 

import ConfigParser
import cmd_parse as cmdp

def set_parse(infile):
    p = ConfigParser.RawConfigParser()
###################### Configuration parameter for Longtrendfilter ###########
    p.add_section('Section GeneralParams')
    p.set('Section GeneralParams','infile','')
    p.set('Section GeneralParams','outfile','')
    p.set('Section GeneralParams','inpath','')
    p.set('Section GeneralParams','outpath','')
    p.set('Section GeneralParams','coljd','')
    p.set('Section GeneralParams','colmag','')
    p.set('Section GeneralParams','lcflag','0')
    p.set('Section GeneralParams','Free_params','')
    p.add_section('Section TransitParams')
    p.set('Section TransitParams','period','')
    p.set('Section TransitParams','epoch','')
    p.set('Section TransitParams','rpstar','')
    p.set('Section TransitParams','Tdur','')
    p.set('Section TransitParams','qgress','')
    p.set('Section TransitParams','inc','')
    p.set('Section TransitParams','u1','')
    p.set('Section TransitParams','u2','')
    p.set('Section TransitParams','starf','')
    p.set('Section TransitParams','planetf','')
    p.set('Section TransitParams','alpha','')
    p.set('Section TransitParams','stdmag','')
    p.add_section('Section FitParams')
    p.set('Section FitParams','period','')
    p.set('Section FitParams','period_err','')
    p.set('Section FitParams','t0','')
    p.set('Section FitParams','t0_err','')
    p.set('Section FitParams','rp','')
    p.set('Section FitParams','rp_err','')
    p.set('Section FitParams','inc','')
    p.set('Section FitParams','inc_err','')
    p.set('Section FitParams','u1','')
    p.set('Section FitParams','u1_err','')
    p.set('Section FitParams','u2','')
    p.set('Section FitParams','u2_err','')
    p.set('Section FitParams','starf','')
    p.set('Section FitParams','starf_err','')
    p.set('Section FitParams','planetf','')
    p.set('Section FitParams','planetf_err','')
    p.set('Section FitParams','rsum','')
    p.set('Section FitParams','rsum_err','')


#################### The General Setting part #############################
	
    p.add_section('This is a configure file for xxx package')

    with open (infile,'wb') as configfile:
        p.write(configfile)	
        return 

def tran_parse(infile,name):

    p = ConfigParser.RawConfigParser()
    try:
        p.read(infile)
        var=p.get('Section TransitParams',name)

    except IOError:
        print '\n'
    except ConfigParser.MissingSectionHeaderError:
        print 'Error: the section TransitParams is missing, excute set_parse to see example.cfg\n'
    except ConfigParser.NoOptionError:
        print 'Error: the option %s is missing, excute set_parse to see example.cfg\n' % name
    return var 

def File_parse(infile,name):

    p = ConfigParser.RawConfigParser()
    try:
        p.read(infile)
        var=p.get('Section GeneralParams',name)

    except IOError:
        print '\n'
    except ConfigParser.MissingSectionHeaderError:
        print 'Error: the section GeneralParams is missing, excute set_parse to see example.cfg\n'
    except ConfigParser.NoOptionError:
        print 'Error: the option %s is missing, excute set_parse to see example.cfg\n' % name
    return var 


def Fit_parse(infile,name):

    p = ConfigParser.RawConfigParser()
    try:
        p.read(infile)
        var=p.get('Section FitParams',name)

    except IOError:
        print '\n'
    except ConfigParser.MissingSectionHeaderError:
        print 'Error: the section FitParams is missing, excute set_parse to see example.cfg\n'
    except ConfigParser.NoOptionError:
        print 'Error: the option %s is missing, excute set_parse to see example.cfg\n' % name
    return var 


if __name__=='__main__':
    options=cmdp.runcfg()
    if(options.eflag):	
        infile=options.infile
        set_parse(infile)

