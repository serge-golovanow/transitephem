'''
Ephemeris calculating tool that uses transit data from exoplanets.org
and astrometric calculations by PyEphem to tell you what transits you'll
be able to observe from your observatory in the near future.

Exoplanets.org citation: Wright et al.2011
http://arxiv.org/pdf/1012.5676v3.pdf

Core developer: Brett Morris
Edits : Serge Golovanow

'''
import math #pour math.pi
#from datetime import datetime #unused
from pytz import timezone
import ephem	 ## PyEphem module
import numpy as np
from glob import glob
import time
from astropy.time import Time
import os.path
import pickle
import sys
from urllib.request import urlopen
import re
from timezonefinder import TimezoneFinder

tf = TimezoneFinder()

rootdir = './outputs/'
exodbPath = './outputs/'#rootdir+'exodb'
sortable = '/* https://www.kryogenix.org/code/browser/sorttable/ */ var stIsIE=!1;if(sorttable={init:function(){arguments.callee.done||(arguments.callee.done=!0,_timer&&clearInterval(_timer),document.createElement&&document.getElementsByTagName&&(sorttable.DATE_RE=/^(\\d\\d?)[\\/\\.-](\\d\\d?)[\\/\\.-]((\\d\\d)?\\d\\d)$/,forEach(document.getElementsByTagName("table"),function(t){-1!=t.className.search(/\\bsortable\\b/)&&sorttable.makeSortable(t)})))},makeSortable:function(t){if(0==t.getElementsByTagName("thead").length&&(the=document.createElement("thead"),the.appendChild(t.rows[0]),t.insertBefore(the,t.firstChild)),null==t.tHead&&(t.tHead=t.getElementsByTagName("thead")[0]),1==t.tHead.rows.length){sortbottomrows=[];for(var e=0;e<t.rows.length;e++)-1!=t.rows[e].className.search(/\\bsortbottom\\b/)&&(sortbottomrows[sortbottomrows.length]=t.rows[e]);if(sortbottomrows){null==t.tFoot&&(tfo=document.createElement("tfoot"),t.appendChild(tfo));for(e=0;e<sortbottomrows.length;e++)tfo.appendChild(sortbottomrows[e]);delete sortbottomrows}headrow=t.tHead.rows[0].cells;for(e=0;e<headrow.length;e++)headrow[e].className.match(/\\bsorttable_nosort\\b/)||(mtch=headrow[e].className.match(/\\bsorttable_([a-z0-9]+)\\b/),mtch&&(override=mtch[1]),mtch&&"function"==typeof sorttable["sort_"+override]?headrow[e].sorttable_sortfunction=sorttable["sort_"+override]:headrow[e].sorttable_sortfunction=sorttable.guessType(t,e),headrow[e].sorttable_columnindex=e,headrow[e].sorttable_tbody=t.tBodies[0],dean_addEvent(headrow[e],"click",sorttable.innerSortFunction=function(t){if(-1!=this.className.search(/\\bsorttable_sorted\\b/))return sorttable.reverse(this.sorttable_tbody),this.className=this.className.replace("sorttable_sorted","sorttable_sorted_reverse"),this.removeChild(document.getElementById("sorttable_sortfwdind")),sortrevind=document.createElement("span"),sortrevind.id="sorttable_sortrevind",sortrevind.innerHTML=stIsIE?\'&nbsp<font face="webdings">5</font>\':"&nbsp;&#x25B4;",void this.appendChild(sortrevind);if(-1!=this.className.search(/\\bsorttable_sorted_reverse\\b/))return sorttable.reverse(this.sorttable_tbody),this.className=this.className.replace("sorttable_sorted_reverse","sorttable_sorted"),this.removeChild(document.getElementById("sorttable_sortrevind")),sortfwdind=document.createElement("span"),sortfwdind.id="sorttable_sortfwdind",sortfwdind.innerHTML=stIsIE?\'&nbsp<font face="webdings">6</font>\':"&nbsp;&#x25BE;",void this.appendChild(sortfwdind);theadrow=this.parentNode,forEach(theadrow.childNodes,function(t){1==t.nodeType&&(t.className=t.className.replace("sorttable_sorted_reverse",""),t.className=t.className.replace("sorttable_sorted",""))}),sortfwdind=document.getElementById("sorttable_sortfwdind"),sortfwdind&&sortfwdind.parentNode.removeChild(sortfwdind),sortrevind=document.getElementById("sorttable_sortrevind"),sortrevind&&sortrevind.parentNode.removeChild(sortrevind),this.className+=" sorttable_sorted",sortfwdind=document.createElement("span"),sortfwdind.id="sorttable_sortfwdind",sortfwdind.innerHTML=stIsIE?\'&nbsp<font face="webdings">6</font>\':"&nbsp;&#x25BE;",this.appendChild(sortfwdind),row_array=[],col=this.sorttable_columnindex,rows=this.sorttable_tbody.rows;for(var e=0;e<rows.length;e++)row_array[row_array.length]=[sorttable.getInnerText(rows[e].cells[col]),rows[e]];row_array.sort(this.sorttable_sortfunction),tb=this.sorttable_tbody;for(e=0;e<row_array.length;e++)tb.appendChild(row_array[e][1]);delete row_array}))}},guessType:function(t,e){sortfn=sorttable.sort_alpha;for(var r=0;r<t.tBodies[0].rows.length;r++)if(text=sorttable.getInnerText(t.tBodies[0].rows[r].cells[e]),""!=text){if(text.match(/^-?[£$¤]?[\\d,.]+%?$/))return sorttable.sort_numeric;if(possdate=text.match(sorttable.DATE_RE),possdate){if(first=parseInt(possdate[1]),second=parseInt(possdate[2]),first>12)return sorttable.sort_ddmm;if(second>12)return sorttable.sort_mmdd;sortfn=sorttable.sort_ddmm}}return sortfn},getInnerText:function(t){if(!t)return"";if(hasInputs="function"==typeof t.getElementsByTagName&&t.getElementsByTagName("input").length,null!=t.getAttribute("sorttable_customkey"))return t.getAttribute("sorttable_customkey");if(void 0!==t.textContent&&!hasInputs)return t.textContent.replace(/^\\s+|\\s+$/g,"");if(void 0!==t.innerText&&!hasInputs)return t.innerText.replace(/^\\s+|\\s+$/g,"");if(void 0!==t.text&&!hasInputs)return t.text.replace(/^\\s+|\\s+$/g,"");switch(t.nodeType){case 3:if("input"==t.nodeName.toLowerCase())return t.value.replace(/^\\s+|\\s+$/g,"");case 4:return t.nodeValue.replace(/^\\s+|\\s+$/g,"");case 1:case 11:for(var e="",r=0;r<t.childNodes.length;r++)e+=sorttable.getInnerText(t.childNodes[r]);return e.replace(/^\\s+|\\s+$/g,"");default:return""}},reverse:function(t){newrows=[];for(var e=0;e<t.rows.length;e++)newrows[newrows.length]=t.rows[e];for(e=newrows.length-1;e>=0;e--)t.appendChild(newrows[e]);delete newrows},sort_numeric:function(t,e){return aa=parseFloat(t[0].replace(/[^0-9.-]/g,"")),isNaN(aa)&&(aa=0),bb=parseFloat(e[0].replace(/[^0-9.-]/g,"")),isNaN(bb)&&(bb=0),aa-bb},sort_alpha:function(t,e){return t[0]==e[0]?0:t[0]<e[0]?-1:1},sort_ddmm:function(t,e){return mtch=t[0].match(sorttable.DATE_RE),y=mtch[3],m=mtch[2],d=mtch[1],1==m.length&&(m="0"+m),1==d.length&&(d="0"+d),dt1=y+m+d,mtch=e[0].match(sorttable.DATE_RE),y=mtch[3],m=mtch[2],d=mtch[1],1==m.length&&(m="0"+m),1==d.length&&(d="0"+d),dt2=y+m+d,dt1==dt2?0:dt1<dt2?-1:1},sort_mmdd:function(t,e){return mtch=t[0].match(sorttable.DATE_RE),y=mtch[3],d=mtch[2],m=mtch[1],1==m.length&&(m="0"+m),1==d.length&&(d="0"+d),dt1=y+m+d,mtch=e[0].match(sorttable.DATE_RE),y=mtch[3],d=mtch[2],m=mtch[1],1==m.length&&(m="0"+m),1==d.length&&(d="0"+d),dt2=y+m+d,dt1==dt2?0:dt1<dt2?-1:1},shaker_sort:function(t,e){for(var r=0,o=t.length-1,n=!0;n;){n=!1;for(var s=r;s<o;++s)if(e(t[s],t[s+1])>0){var a=t[s];t[s]=t[s+1],t[s+1]=a,n=!0}if(o--,!n)break;for(s=o;s>r;--s)if(e(t[s],t[s-1])<0){a=t[s];t[s]=t[s-1],t[s-1]=a,n=!0}r++}}},document.addEventListener&&document.addEventListener("DOMContentLoaded",sorttable.init,!1),/WebKit/i.test(navigator.userAgent))var _timer=setInterval(function(){/loaded|complete/.test(document.readyState)&&sorttable.init()},10);function dean_addEvent(t,e,r){if(t.addEventListener)t.addEventListener(e,r,!1);else{r.$$guid||(r.$$guid=dean_addEvent.guid++),t.events||(t.events={});var o=t.events[e];o||(o=t.events[e]={},t["on"+e]&&(o[0]=t["on"+e])),o[r.$$guid]=r,t["on"+e]=handleEvent}}function removeEvent(t,e,r){t.removeEventListener?t.removeEventListener(e,r,!1):t.events&&t.events[e]&&delete t.events[e][r.$$guid]}function handleEvent(t){var e=!0;t=t||fixEvent(((this.ownerDocument||this.document||this).parentWindow||window).event);var r=this.events[t.type];for(var o in r)this.$$handleEvent=r[o],!1===this.$$handleEvent(t)&&(e=!1);return e}function fixEvent(t){return t.preventDefault=fixEvent.preventDefault,t.stopPropagation=fixEvent.stopPropagation,t}window.onload=sorttable.init,dean_addEvent.guid=1,fixEvent.preventDefault=function(){this.returnValue=!1},fixEvent.stopPropagation=function(){this.cancelBubble=!0},Array.forEach||(Array.forEach=function(t,e,r){for(var o=0;o<t.length;o++)e.call(r,t[o],o,t)}),Function.prototype.forEach=function(t,e,r){for(var o in t)void 0===this.prototype[o]&&e.call(r,t[o],o,t)},String.forEach=function(t,e,r){Array.forEach(t.split(""),function(o,n){e.call(r,o,n,t)})};var forEach=function(t,e,r){if(t){var o=Object;if(t instanceof Function)o=Function;else{if(t.forEach instanceof Function)return void t.forEach(e,r);"string"==typeof t?o=String:"number"==typeof t.length&&(o=Array)}o.forEach(t,e,r)}};'
css = '#eph td,.headinfo,h1{text-align:center}body{font-family:helvetica,arial,sans-serif;padding:0em;background-color:grey;margin-left:auto;margin-right:auto;width:65em;margin-bottom:auto}.headinfo{font-size:80%;margin-left:0}table.daynight,table.sortable{margin-left:auto;margin-right:auto}.small{font-size:8pt;font-style:normal}.linkList{list-style-type:none;padding-left:4em}#eph td{padding-left:.8em;padding-right:.8em}#eph tr:hover{background-color:#D8D8D8}#eph th{background-color:#386890;color:#FFF}#eph tbody tr:nth-child(odd){background-color:#EAF5F5}#eph tr:nth-child(even){background-color:#FFF}#eph tbody tr:nth-child(odd):hover,#eph tr:nth-child(even):hover{background-color:#BDFCC9}table.sortable{padding:1em;font-size:12pt}table.daynight{font-size:10pt;background-color:#EAF5F5}div{display:block;padding:1em;border:solid silver;border-radius:.5em;background-color:#fff}a{color:#06C;text-decoration:none;font-size:105%}a:hover{color:#009;text-decoration:underline}p{margin-left:2em;margin-right:2em}li{line-height:1.2em} '




def gd2jd(gdlist):
    '''
    Parameters
    ----------
    gdlist : list
        Gregorian date in list format - [year,month,day,hour,minute,second]
    Returns
    -------
    jd : float
        Julian date
    '''
    
    # Convert the input list into a string in "ISO" format for astropy
    gdstring = "%i-%i-%i %s:%s:%s" % tuple(gdlist)
    return Time(gdstring, format='iso', scale='utc').jd

def jd2gd(jd):
    '''
    Parameters
    ----------
    jd : float
        Time in julian date
    Returns
    -------
    gdlist : list
        Gregorian date in list format - [year,month,day,hour,minute,second]
    '''
    string = Time(jd, format='jd', scale='utc').iso
    return tuple(map(int,re.findall(r"[\w']+",string)[:-1]))


def downloadAndPickle():
    pklDatabaseName = os.path.join(exodbPath,'exoplanetDB.pkl')	 ## Name of exoplanet database C-pickle
    pklDatabasePaths = glob(pklDatabaseName)   ## list of files with the name pklDatabaseName in cwd
    csvDatabaseName = os.path.join(exodbPath,'exoplanets.csv')  ## Path to the text file saved from exoplanets.org
    csvDatabasePaths = glob(csvDatabaseName)

    '''First, check if there is an internet connection.'''

    '''If there's a previously archived database pickle in this current working 
        directory then use it, if not, grab the data from exoplanets.org in one big CSV file and make one.
        If the old archive is >14 days old, grab a fresh version of the database from exoplanets.org.
        '''
    if csvDatabasePaths == []:
        #print('No local copy of exoplanets.org database. Downloading one...')
        rawCSV = urlopen('http://www.exoplanets.org/csv-files/exoplanets.csv').read()
        saveCSV = open(csvDatabaseName,'w')
        saveCSV.write(str(rawCSV))
        saveCSV.close()
    else: 
        '''If the local copy of the exoplanets.org database is >14 days old, download a new one'''
        secondsSinceLastModification = time.time() - os.path.getmtime(csvDatabaseName) ## in seconds
        daysSinceLastModification = secondsSinceLastModification/(60*60*24*30)
        if daysSinceLastModification > 7:
            #print('Your local copy of the exoplanets.org database is >14 days old. Downloading a fresh one...')
            rawCSV = urlopen('http://www.exoplanets.org/csv-files/exoplanets.csv').read()
            saveCSV = open(csvDatabaseName,'w')
            saveCSV.write(rawCSV)
            saveCSV.close()
        #else: #print("Your local copy of the exoplanets.org database is <14 days old. That'll do.")

    if len(pklDatabasePaths) == 0:
        print('Parsing '+os.path.split(csvDatabaseName)[1]+', the CSV database from exoplanets.org...')
        rawTable = open(csvDatabaseName).read().splitlines()
        labels = rawTable[0].split(',')
        #labelUnits = rawTable[1].split(',')
        #rawTableArray = np.zeros([len(rawTable),len(labels)])
        exoplanetDB = {}
        planetNameColumn = np.arange(len(labels))[np.array(labels,dtype=str)=='NAME'][0]
        for row in range(1,len(rawTable)): 
            splitRow = rawTable[row].split(',')
            exoplanetDB[splitRow[planetNameColumn]] = {}	## Create dictionary for this row's planet
            for col in range(0,len(splitRow)):
                exoplanetDB[splitRow[planetNameColumn]][labels[col]] = splitRow[col]
        
        output = open(pklDatabaseName,'wb')
        pickle.dump(exoplanetDB,output)
        output.close()
    else: 
        #print('Using previously parsed database from exoplanets.org...')
        inputFile = open(pklDatabaseName,'rb')
        exoplanetDB = pickle.load(inputFile)
        inputFile.close()
    
    return exoplanetDB

def calculateEphemerides(parFile,conf):
    '''
        :INPUTS:
        parFile	 --	  path to the parameter file
        '''


    '''Parse the observatory .par file'''
    parFileText = open(os.path.join(os.path.dirname(__file__),parFile),'r').read().splitlines()

    def returnBool(value):
        '''Return booleans from strings'''
        if value.upper().strip() == 'TRUE': return True
        elif value.upper().strip() == 'FALSE': return False
    if hasattr(sys, 'real_prefix'):
        show_lt = float(0)
    for line in parFileText:
        parameter = line.split('=')[0]
        if len(line.split('=')) > 1:
            value = line.split('=')[1].strip()
            if parameter == 'name': observatory_name = value
            elif parameter == 'latitude': observatory_latitude = value
            elif parameter == 'longitude': observatory_longitude = value
            elif parameter == 'elevation': observatory_elevation = float(value)
#            elif parameter == 'temperature': observatory_temperature = float(value)
            elif parameter == 'min_horizon': observatory_minHorizon = value
            elif parameter == 'start_date': startSem = gd2jd(eval(value))
            elif parameter == 'end_date': endSem = gd2jd(eval(value))
            elif parameter == 'mag_limit': mag_limit = float(value)
            elif parameter == 'band': band = value
            elif parameter == 'depth_limit': depth_limit = float(value)
            elif parameter == 'calc_transits': calcTransits = returnBool(value)
            elif parameter == 'calc_eclipses': calcEclipses = returnBool(value)
            elif parameter == 'html_out': htmlOut = returnBool(value)
#            elif parameter == 'text_out': textOut = returnBool(value)
            elif parameter == 'twilight': twilightType = value
            elif parameter == 'max_duration': max_duration = int(value)
            elif parameter == 'show_lt': show_lt = float(value)
    #from oscaar.extras.knownSystemParameters import getLatestParams
    exoplanetDB = downloadAndPickle()
    
    observatory_latitude = conf['latitude']
    observatory_longitude = conf['longitude']
    observatory_elevation = float(conf['elevation'])
    observatory_minHorizon = conf['horizon']
    startSem = gd2jd(eval(conf['start'])+ (0,0,0))
    endSem = gd2jd(eval(conf['end'])+ (0,0,0))
    if endSem-startSem>31: endSem = startSem+31 #max 31jrs calculés
    mag_limit = float(conf['mag'])
    band = 'V'
    depth_limit = float(conf['depth'])
    twilightType = conf['twilight']
    max_duration = int(conf['duration'])
    
    htmlOut = True
    calcTransits = True
    calcEclipses = False
    show_lt=0

    ''' Set up observatory parameters '''
    #ephem.degrees
    
    observatory = ephem.Observer()
    observatory.lat =  observatory_latitude#'38:58:50.16'	## Input format-  deg:min:sec  (type=str)
    observatory.long = observatory_longitude#'-76:56:13.92' ## Input format-  deg:min:sec  (type=str)
    observatory.elevation = observatory_elevation   # m
    #observatory.temp = observatory_temperature	  ## Celsius 
    observatory.horizon = observatory_minHorizon	## Input format-  deg:min:sec  (type=str)
    
    tzlat = observatory.lat*180/math.pi
    tzlng = observatory.long*180/math.pi
    tz = tf.timezone_at(lng=float(tzlng), lat=tzlat)

    def trunc(f, n):
        '''Truncates a float f to n decimal places without rounding'''
        slen = len('%.*f' % (n, f))
        return str(f)[:slen]

    def RA(planet):
        '''Type: str, Units:  hours:min:sec'''
        return exoplanetDB[planet]['RA_STRING']
    def dec(planet):
        '''Type: str, Units:  deg:min:sec'''
        return exoplanetDB[planet]['DEC_STRING']
    def period(planet):
        '''Units:  days'''
        return np.float64(exoplanetDB[planet]['PER'])
    def epoch(planet):
        '''Tc at mid-transit. Units:  days'''
        if exoplanetDB[planet]['TT'] == '': return 0.0
        else: return np.float64(exoplanetDB[planet]['TT'])
    def duration(planet):
        '''Transit/eclipse duration. Units:  days'''
        if exoplanetDB[planet]['T14'] == '': return 0.0
        else: return float(exoplanetDB[planet]['T14'])
    def V(planet):
        '''V mag'''
        if exoplanetDB[planet]['V'] == '': return 0.0
        else: return float(exoplanetDB[planet]['V'])
    def KS(planet):
        '''KS mag'''
        if exoplanetDB[planet]['KS'] == '': return 0.0
        else: return float(exoplanetDB[planet]['KS'])
    
    def bandMagnitude(planet):
        if band.upper() == 'V':
            return V(planet)
        elif band.upper() == 'K':
            return KS(planet)
    def depth(planet):
        '''Transit depth'''
        if exoplanetDB[planet]['DEPTH'] == '': return 0.0
        else: return float(exoplanetDB[planet]['DEPTH'])

    def transitBool(planet):
        '''True if exoplanet is transiting, False if detected by other means'''
        if exoplanetDB[planet]['TRANSIT'] == '0': return 0
        elif exoplanetDB[planet]['TRANSIT'] == '1': return 1
    ########################################################################################
    ########################################################################################

    def datestr2list(datestr):
        ''' Take strings of the form: "2013/1/18 20:08:18" and return them as a
            tuple of the same parameters'''
        year,month,others = datestr.split('/')
        day, time = others.split(' ')
        hour,minute,sec = time.split(':')
        return (int(year),int(month),int(day),int(hour),int(minute),int(sec))

    def list2datestr(inList):
        '''Converse function to datestr2list'''
        inList = list(map(str,inList))
        return inList[0]+'/'+inList[1]+'/'+inList[2]+' '+inList[3].zfill(2)+':'+inList[4].zfill(2)+':'+inList[5].zfill(2)

    def list2datestrCSV(inList):
        '''Converse function to datestr2list'''
        inList = list(map(str,inList))
        #print inList
        return inList[0]+'/'+inList[1]+'/'+inList[2]+','+inList[3].zfill(2)+':'+inList[4].zfill(2)+':'+inList[5].zfill(2)


    def list2datestrHTML(inList,alt,direction):
        '''Converse function to datestr2list'''
        time = gd2jd(list(inList))
        etime=ephem.Date(inList) 
        #tz=pytz.timezone('Europe/Paris')
        dt = etime.datetime().replace(tzinfo=timezone('UTC'))
        if tz != None:
            dt = dt.astimezone(timezone(tz))
            tl = ephem.localtime(etime)
            #tl = tl.astimezone(pytz.timezone('Europe/Berlin'))
            tl = dt.time()
            tl = str(tl.hour).zfill(2) + ':' + str(tl.minute).zfill(2)
        else:
            tl = 'inconnue'
        inList = list(map(str,inList))
        altsun = str(trunc(altSoleil(time),0))
        #return inList[1].zfill(2)+'/'+inList[2].zfill(2)+'<br />'+inList[3].zfill(2)+':'+inList[4].zfill(2)
        return inList[2].zfill(2)+'/'+inList[1].zfill(2)+', <abbr title="Heure l&eacute;gale '+str(tl)+',&#10;hauteur Soleil '+altsun+'&#xB0;">'+inList[3].zfill(2)+':'+inList[4].split('.')[0].zfill(2)+'</abbr><br /> '+alt+'&deg; '+direction

    def list2datestrHTML_UTnoaltdir(inList,alt,direction):
        '''Converse function to datestr2list'''
        inList = list(map(str,inList))
        #return inList[1].zfill(2)+'/'+inList[2].zfill(2)+'<br />'+inList[3].zfill(2)+':'+inList[4].zfill(2)
        return inList[1].zfill(2)+'/<strong>'+inList[2].zfill(2)+'</strong>, '+inList[3].zfill(2)+':'+inList[4].split('.')[0].zfill(2)

    def list2datestrHTML_LT(inList,alt,direction):
        '''Converse function to datestr2list for daylight savings time'''
        #print "original",inList
        tempDate = ephem.Date(inList)
        inList = ephem.Date(ephem.localtime(tempDate)).tuple()
        #print "converted",lt_inList,'\n'
        inList = list(map(str,inList))
        #return inList[1].zfill(2)+'/'+inList[2].zfill(2)+'<br />'+inList[3].zfill(2)+':'+inList[4].zfill(2)
        return inList[1].zfill(2)+'/<strong>'+inList[2].zfill(2)+'</strong>, '+inList[3].zfill(2)+':'+inList[4].split('.')[0].zfill(2)+'<br /> '+alt+'&deg; '+direction

    def simbadURL(planet):
        if exoplanetDB[planet]['SIMBADURL'] == '': return 'http://simbad.harvard.edu/simbad/'
        else: return exoplanetDB[planet]['SIMBADURL']

    def RADecHTML(planet):
        return '<a href="'+simbadURL(planet)+'">'+RA(planet).split('.')[0]+'<br />'+dec(planet).split('.')[0]+'</a>'

    def constellation(planet):
        return exoplanetDB[planet]['Constellation']
    
    def url(planet):
        return 'http://www.exoplanets.org/detail/'+planet
    
    def etd(planet):
        return '<abbr title="Exoplanet Transit Database"><a href="http://var2.astro.cz/ETD/etd.php?STARNAME='+exoplanetDB[planet]['STAR']+'&PLANET='+exoplanetDB[planet]['COMP']+'">ETD</a></abbr>'
    
    def epe(planet):
        return '<abbr title="Extrasolar Planets Encyclopaedia"><a href="http://exoplanet.eu/catalog/'+planet+'">EPE</a></abbr>'
    
    def oec(planet):
        return '<abbr title="Open Exoplanet Catalog"><a href="http://www.openexoplanetcatalogue.com/planet/'+planet+'">OEC</a></abbr>'
    
    def nasa(planet):
        return '<a href="https://exoplanetarchive.ipac.caltech.edu/cgi-bin/DisplayOverview/nph-DisplayOverview?objname='+exoplanetDB[planet]['STAR']+'+'+exoplanetDB[planet]['COMP']+'">NASA</a>'

    def orbitReference(planet):
        return exoplanetDB[planet]['TRANSITURL']

    def orbitReferenceYear(planet):
        '''ORBREF returns the citation in the format "<first author> <year>", so parse and return just the year'''
        return exoplanetDB[planet]['ORBREF'].split()[1]

    def nameWithLink(planet):
        return '<a href="'+url(planet)+'"><b>'+planet+'</b></a>'+'<span class="small"><br />'+etd(planet)+', '+oec(planet)+', '+epe(planet)+', '+nasa(planet)+'</span>'

    def mass(planet):
        if exoplanetDB[planet]['MASS'] == '': return '---'
        else: return trunc(float(exoplanetDB[planet]['MASS']),2)

    def semimajorAxis(planet):
        #return trunc(0.004649*float(exoplanetDB[planet]['AR'])*float(exoplanetDB[planet]['RSTAR']),3)   ## Convert from solar radii to AU
        return trunc(float(exoplanetDB[planet]['SEP']),3)

    def radius(planet):
        if exoplanetDB[planet]['R'] == '': return '---'
        else: return trunc(float(exoplanetDB[planet]['R']),2) ## Convert from solar radii to Jupiter radii

    def midTransit(Tc, P, start, end):
        '''Calculate mid-transits between Julian Dates start and end, using a 2500 
            orbital phase kernel since T_c (for 2 day period, 2500 phases is 14 years)
            '''
        Nepochs = np.arange(0,2500,dtype=np.float64)
        transitTimes = Tc + P*Nepochs
        transitTimesInSem = transitTimes[(transitTimes < end)*(transitTimes > start)]
        return transitTimesInSem

    def midEclipse(Tc, P, start, end):
        '''Calculate mid-eclipses between Julian Dates start and end, using a 2500 
            orbital phase kernel since T_c (for 2 day period, 2500 phases is 14 years)
            '''
        Nepochs = np.arange(0,2500,dtype=np.float64)
        transitTimes = Tc + P*(0.5 + Nepochs)
        transitTimesInSem = transitTimes[(transitTimes < end)*(transitTimes > start)]
        return transitTimesInSem

    '''Choose which planets from the database to include in the search, 
        assemble a list of them.'''
    planets = []
    for planet in exoplanetDB:
        if bandMagnitude(planet) != 0.0 and depth(planet) != 0.0 and float(bandMagnitude(planet)) <= mag_limit and \
           float(depth(planet)) >= depth_limit and duration(planet) != 0.0 and transitBool(planet) and duration(planet)*24 <= max_duration:
            planets.append(planet)

    if calcTransits: transits = {}
    if calcEclipses: eclipses = {}
    for day in np.arange(startSem,endSem+1):
        if calcTransits: transits[str(day)] = []
        if calcEclipses: eclipses[str(day)] = []
    planetsNeverUp = []


    def azToDirection(az):
        az = float(az)
        if (az >= 0 and az < 22.5) or (az >= 337.5 and az < 360): lettre='N'
        elif az >= 22.5 and az < 67.5:  lettre='NE'
        elif az >= 67.5 and az < 112.5:  lettre='E'
        elif az >= 112.5 and az < 157.5:  lettre='SE'
        elif az >= 157.5 and az < 202.5:  lettre='S'
        elif az >= 202.5 and az < 247.5:  lettre='SW'
        elif az >= 247.5 and az < 292.5:  lettre='W'	
        elif az >= 292.5 and az < 337.5:  lettre='NW'
        return '<abbr title="azimuth '+str(trunc(az,0))+'&#xB0;">'+lettre+'</abbr>'

    def ingressEgressAltAz(planet,observatory,ingress,egress,mid):
        altitudes = []
        directions = []
        for time in [ingress,egress,mid]:
            observatory.date = list2datestr(jd2gd(time))
            star = ephem.FixedBody()
            star._ra = ephem.hours(RA(planet))
            star._dec = ephem.degrees(dec(planet))
            star.compute(observatory)
            altitudes.append(str(ephem.degrees(star.alt)).split(":")[0])
            directions.append(azToDirection(str(ephem.degrees(star.az)).split(":")[0]))
        ingressAlt,egressAlt,midAlt = altitudes
        ingressDir,egressDir,midDir = directions
        return ingressAlt,ingressDir,egressAlt,egressDir,midAlt,midDir

    def aboveHorizonForEvent(planet,observatory,ingress,egress):
        altitudes = []
        for time in [ingress,egress]:
            observatory.date = list2datestr(jd2gd(time))
            star = ephem.FixedBody()
            star._ra = ephem.hours(RA(planet))
            star._dec = ephem.degrees(dec(planet))
            star.compute(observatory)
            #altitudes.append(str(ephem.degrees(star.alt)).split(":")[0])
            altitudes.append(float(repr(star.alt))/(2*np.pi) * 360)	## Convert altitudes to degrees
        #if altitudes[0] > 0 and altitudes[1] > 0: return True
        if altitudes[0] > float(ephem.degrees(observatory_minHorizon))*(180/np.pi) and altitudes[1] > float(ephem.degrees(observatory_minHorizon))*(180/np.pi): return True
        else: return False
        
    def lunaison(time,planet):
        observatory.date = list2datestr(jd2gd(time))
        moon = ephem.Moon()
        moon.compute(observatory)
        lunes = ['&#127761;','&#127762;','&#127763;','&#127764;','&#127765;','&#127766;','&#127767;','&#127768;']
        nnm = ephem.next_new_moon(observatory.date)  
        pnm = ephem.previous_new_moon(observatory.date)  
        lun = (observatory.date - pnm) / (nnm - pnm)     
        if lun <=0.0625 or lun>0.9375: ilun=0
        elif lun <=0.1875: ilun=1
        elif lun <=0.3125: ilun=2
        elif lun <=0.4375: ilun=3        
        elif lun <=0.5625: ilun=4        
        elif lun <=0.6875: ilun=5        
        elif lun <=0.8125: ilun=6        
        elif lun <=0.9375: ilun=7        
        #elif lun <=0.3125: ilun=8  
        #radec= (RA(planet).split('.')[0],dec(planet).split('.')[0])
        star = ephem.readdb("star,f|M|F7,"+RA(planet)+","+dec(planet)+",2.02,2000")
        star.compute();
        
        sepa = str(ephem.separation( moon,  star )).split(':')[0]

        return '<abbr title="&Agrave; mi-transit :&#10;illumination '+trunc(100*moon.moon_phase,0)+'%&#10;&eacute;longation '+str(sepa)+'&#xB0;">'+lunes[ilun]+'</abbr>'
    
    def luneillum(time):
        observatory.date = list2datestr(jd2gd(time))
        moon = ephem.Moon()
        moon.compute(observatory)
        return moon.phase
        
    def altSoleil(time):
        observatory.date = list2datestr(jd2gd(time))
        sun = ephem.Sun()
        sun.compute(observatory)
        return float(repr(sun.alt))/(2*np.pi) * 360

    def eventAfterTwilight(planet,observatory,ingress,egress,twilightType):
        altitudes = []
        for time in [ingress,egress]:
            observatory.date = list2datestr(jd2gd(time))
            sun = ephem.Sun()
            sun.compute(observatory)
            altitudes.append(float(repr(sun.alt))/(2*np.pi) * 360)	## Convert altitudes to degrees
        if altitudes[0] < float(twilightType) and altitudes[1] < float(twilightType): return True
        else: return False

    constellations = {"And":"Androm&eagrave;de","Ant":"La Machine pneumatique","Aps":"L'Oiseau de paradis","Aql":"L'Aigle","Aqr":"Le Verseau","Ara":"L'Autel","Ari":"Le B&eacute;lier","Aur":"Le Cocher","Boo":"Le Bouvier","Cae":"Le Burin","Cam":"La Girafe","Cap":"Le Capricorne","Car":"La Car&eagrave;ne","Cas":"Cassiop&eacute;e","Cen":"Le Centaure","Cep":"C&eacute;ph&eacute;e","Cet":"La Baleine","Cha":"Le Cam&eacute;l&eacute;on","Cir":"Le Compas","CMa":"Le Grand Chien","CMi":"Le Petit Chien","Cnc":"Le Cancer","Col":"La Colombe","Com":"La Chevelure de B&eacute;r&eacute;nice","CrA":"La Couronne australe","CrB":"La Couronne bor&eacute;ale","Crt":"La Coupe","Cru":"La Croix du Sud","Crv":"Le Corbeau","CVn":"Les Chiens de chasse","Cyg":"Le Cygne","Del":"Le Dauphin","Dor":"La Dorade","Dra":"Le Dragon","Equ":"Le Petit Cheval","Eri":"L'Éridan","For":"Le Fourneau","Gem":"Les G&eacute;meaux","Gru":"La Grue","Her":"Hercule","Hor":"L'Horloge","Hya":"L'Hydre","Hyi":"L'Hydre mâle","Ind":"L'Indien","Lac":"Le L&eacute;zard","Leo":"Le Lion","Lep":"Le Li&eagrave;vre","Lib":"La Balance","LMi":"Le Petit Lion","Lup":"Le Loup","Lyn":"Le Lynx","Lyr":"La Lyre","Men":"La Table","Mic":"Le Microscope","Mon":"La Licorne","Mus":"La Mouche","Nor":"La R&eagrave;gle","Oct":"L'Octant","Oph":"Ophiuchus","Ori":"Orion","Pav":"Le Paon","Peg":"P&eacute;gase","Per":"Pers&eacute;e","Phe":"Le Ph&eacute;nix","Pic":"Le Peintre","PsA":"Le Poisson austral","Psc":"Les Poissons","Pup":"La Poupe","Pyx":"La Boussole","Ret":"Le R&eacute;ticule","Scl":"Le Sculpteur","Sco":"Le Scorpion","Sct":"L'Écu de Sobieski","Ser":"Le Serpent","Sex":"Le Sextant","Sge":"La Fl&eagrave;che","Sgr":"Le Sagittaire","Tau":"Le Taureau","Tel":"Le T&eacute;lescope","TrA":"Le Triangle austral","Tri":"Le Triangle","Tuc":"Le Toucan","UMa":"La Grande Ourse","UMi":"La Petite Ourse","Vel":"Les Voiles","Vir":"La Vierge","Vol":"Le Poisson volant","Vul":"Le Petit Renard"}

    for planet in planets:		
        '''Compute all of the coming transits and eclipses for a long time out'''
        allTransitEpochs = midTransit(epoch(planet),period(planet),startSem,endSem)
        allEclipseEpochs = midEclipse(epoch(planet),period(planet),startSem,endSem)
        for day in np.arange(startSem,endSem+1,1.0):
            try:
                '''For each day, gather the transits and eclipses that happen'''
                transitEpochs = allTransitEpochs[(allTransitEpochs <= day+0.5)*(allTransitEpochs > day-0.5)]
                eclipseEpochs = allEclipseEpochs[(allEclipseEpochs <= day+0.5)*(allEclipseEpochs > day-0.5)]
                if calcTransits and len(transitEpochs) != 0:
                    transitEpoch = transitEpochs[0]
                    ingress = transitEpoch-duration(planet)/2
                    egress = transitEpoch+duration(planet)/2
                    
                    ''' Calculate positions of host stars'''
                    star = ephem.FixedBody()
                    star._ra = ephem.hours(RA(planet))
                    star._dec = ephem.degrees(dec(planet))
                    star.compute(observatory)
                    const = ephem.constellation(star)
                    exoplanetDB[planet]['Constellation'] = '<abbr title="'+const[1]+'&#10;'+constellations.get(str(const[0]),'')+'">'+const[0]+'</abbr>'
                    
                    '''If star is above horizon and sun is below horizon during transit/eclipse:'''		
                    if aboveHorizonForEvent(planet,observatory,ingress,egress) and eventAfterTwilight(planet,observatory,ingress,egress,twilightType):
                        ingressAlt,ingressDir,egressAlt,egressDir,midAlt,midDir = ingressEgressAltAz(planet,observatory,ingress,egress,transitEpoch)
                        transitInfo = [planet,transitEpoch,duration(planet)/2,'transit',ingressAlt,ingressDir,egressAlt,egressDir,midAlt,midDir]
                        transits[str(day)].append(transitInfo)		
                if calcEclipses and len(eclipseEpochs) != 0:
                    eclipseEpoch = eclipseEpochs[0]
                    ingress = eclipseEpoch-duration(planet)/2
                    egress = eclipseEpoch+duration(planet)/2
                    
                    ''' Calculate positions of host stars'''
                    star = ephem.FixedBody()
                    star._ra = ephem.hours(RA(planet))
                    star._dec = ephem.degrees(dec(planet))
                    star.compute(observatory)
                    exoplanetDB[planet]['Constellation'] = ephem.constellation(star)[0]
                    
                    if aboveHorizonForEvent(planet,observatory,ingress,egress) and eventAfterTwilight(planet,observatory,ingress,egress,twilightType):
                        ingressAlt,ingressDir,egressAlt,egressDir = ingressEgressAltAz(planet,observatory,ingress,egress)
                        eclipseInfo = [planet,eclipseEpoch,duration(planet)/2,'eclipse',ingressAlt,ingressDir,egressAlt,egressDir]
                        eclipses[str(day)].append(eclipseInfo)	
            
            except ephem.NeverUpError:
                if str(planet) not in planetsNeverUp:
                    print('Note: planet %s is never above the horizon at this observing location.' % (planet))
                    planetsNeverUp.append(str(planet))

    def removeEmptySets(dictionary):
        '''Remove days where there were no transits/eclipses from the transit/eclipse list dictionary. 
            Can't iterate through the transits dictionary with a for loop because it would change length 
            as keys get deleted, so loop through with while loop until all entries are not empty sets'''
        dayCounter = startSem
        while any(dictionary[day] == [] for day in dictionary):	
            if dictionary[str(dayCounter)] == []:
                del dictionary[str(dayCounter)]
            dayCounter += 1

    if calcTransits: removeEmptySets(transits)
    if calcEclipses: removeEmptySets(eclipses)

    events = {}
    
    def mergeDictionaries(dict):
        total=0
        for key in dict:
            if any(key == eventKey for eventKey in events) == False:	## If key does not exist in events,
                if np.shape(dict[key])[0] == 1:	## If new event is the only one on that night, add only it
                    events[key] = [dict[key][0]]
                    total+=1
                else:			## If there were multiple events that night, add them each
                    events[key] = []
                    for event in dict[key]:
                        events[key].append(event)
                        total+=1
            else:
                if np.shape(dict[key])[0] > 1: ## If there are multiple entries to append,
                    for event in dict[key]:
                        events[key].append(event)
                        total+=1
                else:							## If there is only one to add,
                    events[key].append(dict[key][0])
                    total+=1
    if calcTransits: mergeDictionaries(transits)
    if calcEclipses: mergeDictionaries(eclipses)
    
    total=0
    allKeys = list(events.keys())
    for key in allKeys:
        for planet in events[key]:
            total+=1
    


    if htmlOut: 
        '''Write out a text report with the transits/eclipses. Write out the time of 
            ingress, egress, whether event is transit/eclipse, elapsed in time between
            ingress/egress of the temporally isolated events'''
        #report = open(os.path.join(os.path.abspath(rootdir),'eventReport.html'),'w')
        debut = ephem.Date(jd2gd(startSem)).triple() # list2datestr(jd2gd(startSem)).split(' ')[0]
        debut = str(int(debut[2])).zfill(2)+'/'+str(debut[1]).zfill(2)+'/'+str(debut[0])
        fin = ephem.Date(jd2gd(endSem)).triple() # list2datestr(jd2gd(startSem)).split(' ')[0]
        fin = str(int(fin[2])).zfill(2)+'/'+str(fin[1]).zfill(2)+'/'+str(fin[0])
        
        #allKeys = list(events.keys())
        ## http://www.kryogenix.org/code/browser/sorttable/
        htmlheader = '\n'.join([
#                                '<!doctype html>',\
#                                '<html>',\
#                                '	<head>',\
#                                '		<meta http-equiv="content-type" content="text/html; charset=UTF-8" />',\
#                                '		<title>Transits d\'exoplan&egrave;tes</title>',\
#                                '		<style>'+css+'</style>',\
#                                '	   <script type="text/javascript">'+sortable+'</script>',\
#                                '	   <script type="text/javascript">function premiertri() { var myTH = document.getElementsByTagName("th")[2]; sorttable.innerSortFunction.apply(myTH, []); }  </script>',\
#                                '	</head>',\
#                                '	<body onLoad="premiertri();">',\
                                '		<div id="textDiv">',\
                                #'		<h1>Eph&eacute;m&eacute;rides pour '+observatory_latitude+'&#xB0;N '+observatory_longitude+'&#xB0;E</h1>',\
                                '		<p>Observateur &agrave; '+str(observatory.lat).split('.')[0]+'&#xB0;N '+str(observatory.long).split('.')[0]+'&#xB0;E (<abbr title="timezone">tz</abbr> '+str(tz)+'), transits du '+str(debut)+' au '+str(fin)+', &eacute;toile magnitude &le;'+str(mag_limit)+' (bande '+band+'), profondeur &ge;'+str(depth_limit)+'mag, dur&eacute;e &le;'+str(max_duration)+'h, transit &agrave; plus de '+observatory_minHorizon+'&#xB0; au dessus de l\'horizon, Soleil &agrave; moins de '+twilightType+'&#xB0;',\
                                #'	   Click the column headers to sort. ',\
                                #'		<table class="daynight" id="eph">',\
                                #'		<tr><th colspan=2>Toggle Color Scheme</th></tr>',\
                                #'		<tr><td><a href="#" onclick="changeCSS(\'stylesheetEphem.css\', 0);">Day</a></td><td><a href="#" onclick="changeCSS(\'stylesheetEphemDark.css\', 0);">Night</a></td></tr>',\
                                #'		</table>'
                                '      : <b>'+str(total)+'</b> transits observables.</p>',\
                                ])
        
        if show_lt == 0:
            tableheader = '\n'.join([
                                     '\n		<table class="sortable" id="eph">',\
                                     '		<thead><tr> <th>Plan&egrave;te</th>	  <!--<th>Event<br /><span class="small">[Transit/<br />Eclipse]</span></th>-->	<th>D&eacute;but</th> <th>Milieu</th> <th>Fin</th>'+\
                                     '<th>Dur&eacute;e<br />(h)</th> <th>Magnitude<br />&eacute;toile</th> <th>Profondeur<br />(mmag)</th> <th>RA/Dec<br /></span></th> <th>Const.</th> <!--<th>Mass<br />(M<sub>J</sub>)</th>'+\
                                     '<th>Radius<br />(R<sub>J</sub>)</th> <th>Ref. Year</th>--><th>Lune</th></tr></thead> <tbody>'])
        else:
            tableheader = '\n'.join([
                                     '\n        <table class="sortable" id="eph">',\
                                     '        <tr> <th>Planet<br /><span class="small">[Link: Orbit ref.]</span></th>      <!--<th>Event<br /><span class="small">[Transit/<br />Eclipse]</span></th>-->  <th>Ingress <br /><span class="small">(MM/DD<br />HH:MM (LT), Alt., Dir.)</span></th> <th>Egress <br /><span class="small">(MM/DD<br />HH:MM (LT), Alt., Dir.)</span></th>   '+\
                                     '<th>'+band.upper()+'</th> <th>Depth<br />(mag)</th> <th>Duration<br />(hrs)</th> <th>RA/Dec<br /><span class="small">[Link: Simbad ref.]</span></th> <th>Const.</th> <th>Mass<br />(M<sub>J</sub>)</th>'+\
                                     ' <th>Radius<br />(R<sub>J</sub>)</th> <th>Ref. Year</th> <th>Ingress <br /><span class="small">(MM/DD<br />HH:MM (UT))</span></th> <th>Egress <br /><span class="small">(MM/DD<br />HH:MM, (UT))</span></th></tr>'])
    
        
        tablefooter = '\n'.join([
                                 '\n		</tbody> </table>',\
                                 '		',])
        htmlfooter = '\n'.join([
                                '\n		<p class="headinfo">Les colonnes sont triables, cliquez simplement sur l\'en-t&ecirc;te ! &nbsp; - &nbsp; ',\
                                '      Voir aussi les pr&eacute;visions <a href="http://var2.astro.cz/ETD/predictions.php?delka='+observatory_longitude.split(':')[0]+'&submit=submit&sirka='+observatory_latitude.split(':')[0]+'">ETD</a> ou <a href="https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TransitView/nph-visibletbls?dataset=transits&conf&lat='+observatory_latitude.split(':')[0]+'&lon='+observatory_longitude.split(':')[0]+'">NASA</a>',\
                                '		<br /><a href="http://exocalc.space">Exoplanets transits calculator</a> (Serge Golovanow), <a href="https://github.com/bmorris3/transitephem">transitephem</a> (<a href="http://staff.washington.edu/bmmorris/">Brett Morris</a>), <a href="http://rhodesmill.org/pyephem/">PyEphem</a>, <a href="http://www.exoplanets.org/">exoplanets.org</a>.</p>',\
                                '		</div>',\
#                                '	</body>',\
#                                '</html>'\
                                ])
        #report.write(htmlheader)
        #report.write(tableheader)
        

        print(htmlheader)
        print(tableheader)
        
        allKeys = np.array(allKeys)[np.argsort(allKeys)]
        for key in allKeys:
            def writeHTMLtransit():
                indentation = '		'
                if show_lt != 0: 
                    middle = '</td><td>'.join([nameWithLink(planet[0]),list2datestrHTML_LT(jd2gd(float(planet[1]-planet[2])),planet[4],planet[5]),\
                                               list2datestrHTML_LT(jd2gd(float(planet[1]+planet[2])),planet[6],planet[7]),trunc(bandMagnitude(str(planet[0])),2),\
                                               trunc(1*depth(planet[0]),4),trunc(24.0*duration(planet[0]),2),RADecHTML(planet[0]),constellation(planet[0]),\
                                               mass(planet[0]),radius(planet[0]),orbitReferenceYear(planet[0]),list2datestrHTML_UTnoaltdir(jd2gd(float(planet[1]-planet[2])),planet[4],planet[5]),\
                                               list2datestrHTML_UTnoaltdir(jd2gd(float(planet[1]+planet[2])),planet[6],planet[7])])
                else:
                    middle = nameWithLink(planet[0])+'</td>'+\
                             '<td sorttable_customkey="'+str(planet[1]-planet[2])+'">'+list2datestrHTML(jd2gd(float(planet[1]-planet[2])),planet[4],planet[5])+'</td>'+\
                             '<td sorttable_customkey="'+str(planet[1])+'">'+list2datestrHTML(jd2gd(float(planet[1])),planet[8],planet[9])+'</td>'+\
                             '<td sorttable_customkey="'+str(planet[1]+planet[2])+'">'+ list2datestrHTML(jd2gd(float(planet[1]+planet[2])),planet[6],planet[7])+'</td>'+\
                             '<td>'+'<abbr title="'+trunc(24.0*60*duration(planet[0]),1)+' minutes"><a href="'+orbitReference(planet[0])+'">'+trunc(24.0*duration(planet[0]),2)+'</a></abbr>'+'</td>'+\
                             '<td>'+trunc(bandMagnitude(str(planet[0])),2)+'</td>'+\
                             '<td>'+trunc(1000*depth(planet[0]),1)+'</td>'+\
                             '<td>'+RADecHTML(planet[0])+'</td>'+\
                             '<td>'+constellation(planet[0])+'</td>'+\
                             '<td sorttable_customkey="'+str(luneillum(float(planet[1])))+'">'+lunaison(float(planet[1]),planet[0])
    
                    
    
    
                line = indentation+'<tr><td>'+str(middle)+'</td></tr>\n'
                #report.write(line)
                print(line)
            
            def writeHTMLeclipse():
                indentation = '		'
                if show_lt != 0:
                    middle = '</td><td>'.join([nameWithLink(planet[0]),str(planet[3]),list2datestrHTML_LT(jd2gd(float(planet[1]-planet[2])),planet[4],planet[5]),\
                                               list2datestrHTML_LT(jd2gd(float(planet[1]+planet[2])),planet[6],planet[7]),trunc(bandMagnitude(str(planet[0])),2),\
                                               '---',trunc(24.0*duration(planet[0]),2),RADecHTML(planet[0]),constellation(planet[0]),\
                                               mass(planet[0]),radius(planet[0]),orbitReferenceYear(planet[0]),list2datestrHTML_UTnoaltdir(jd2gd(float(planet[1]-planet[2])),planet[4],planet[5]),\
                                               list2datestrHTML_UTnoaltdir(jd2gd(float(planet[1]+planet[2])),planet[6],planet[7])])
                else: 
                    middle = '</td><td>'.join([nameWithLink(planet[0]),str(planet[3]),list2datestrHTML(jd2gd(float(planet[1]-planet[2])),planet[4],planet[5]),\
                                               list2datestrHTML(jd2gd(float(planet[1]+planet[2])),planet[6],planet[7]),trunc(bandMagnitude(str(planet[0])),2),\
                                               '---',trunc(24.0*duration(planet[0]),2),RADecHTML(planet[0]),constellation(planet[0]),\
                                               mass(planet[0]),radius(planet[0]),orbitReferenceYear(planet[0])])

                line = indentation+'<tr><td>'+middle+'</td></tr>\n'
                #report.write(line)
            
            
            if np.shape(events[key])[0] > 1:
                elapsedTime = []
                
                for i in range(1,len(events[key])):
                    nextPlanet = events[key][1]
                    planet = events[key][0]
                    double = False
                    '''If the other planet's ingress is before this one's egress, then'''
                    if ephem.Date(list2datestr(jd2gd(float(nextPlanet[1]-nextPlanet[2])))) -\
                        ephem.Date(list2datestr(jd2gd(float(planet[1]+planet[2])))) > 0.0:
                            double = True
                            elapsedTime.append(ephem.Date(list2datestr(jd2gd(float(nextPlanet[1]-nextPlanet[2])))) - \
                                               ephem.Date(list2datestr(jd2gd(float(planet[1]+planet[2])))))
                    
                    if ephem.Date(list2datestr(jd2gd(float(planet[1]-planet[2])))) - \
                        ephem.Date(list2datestr(jd2gd(float(nextPlanet[1]+nextPlanet[2])))) > 0.0:
                            '''If the other planet's egress is before this one's ingress, then'''
                            double = True
                            elapsedTime.append(ephem.Date(list2datestr(jd2gd(float(planet[1]-planet[2])))) - \
                                               ephem.Date(list2datestr(jd2gd(float(nextPlanet[1]+nextPlanet[2])))))
                
                for planet in events[key]:
                    if calcTransits and planet[3] == 'transit':
                        writeHTMLtransit()
                    if calcEclipses and planet[3] == 'eclipse':
                        writeHTMLeclipse()		  
            elif np.shape(events[key])[0] == 1:
                planet = events[key][0]
                if calcTransits and planet[3] == 'transit':
                    writeHTMLtransit()
                if calcEclipses and planet[3] == 'eclipse':
                    writeHTMLeclipse()
        #report.write(tablefooter)
        #report.write(htmlfooter)
        #report.close()
        print(tablefooter)
        print(htmlfooter)
    #print exoplanetDB['HD 209458 b']

