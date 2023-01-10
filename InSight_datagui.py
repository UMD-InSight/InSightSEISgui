from tkinter import *
from tkinter import ttk, filedialog
from tkinter.filedialog import askopenfile
from PIL import Image, ImageTk

# Import these for the catalog generation
import csv
import requests
import xml.etree.ElementTree as ET

# Define the name and size of the interface
root = Tk()
root.title('InSight Event Data Processor')
root.geometry("1200x750")

################################################################################
# Put the functions here
################################################################################

############################################
# CATALOG MAKER TO SEARCH FOR THE EVENTS
############################################

def makecatalog():
    global catbox
    global Amode
    global Bmode
    global Cmode
    global Dmode
    global LFmode
    global BBmode
    global HZmode
    global HFmode
    global VFmode

    Avar=Amode.get()
    Bvar=Bmode.get()
    Cvar=Cmode.get()
    Dvar=Dmode.get()

    LFvar=LFmode.get()
    BBvar=BBmode.get()
    HZvar=HZmode.get()
    HFvar=HFmode.get()
    VFvar=VFmode.get()

    catbox.delete(0,END)
    ns = {"q":"http://quakeml.org/xmlns/quakeml/1.2",
    "d":"http://quakeml.org/xmlns/bed/1.2",
    "catalog":"http://anss.org/xmlns/catalog/0.1",
    "tensor":"http://anss.org/xmlns/tensor/0.1",
    "mars":"http://quakeml.org/xmlns/bed/1.2/mars"}
    #
    tree = ET.parse('events_mars_extended_multiorigin_v12_2022-07-01.xml')
    root = tree.getroot()
    eventlist = root.findall('d:eventParameters',ns)
    #
    for ep in eventlist:
            xevents = ep.findall('d:event',ns)
    #
    import os
    if os.path.exists("SeismicCatalog"):
      os.remove("SeismicCatalog")
    #
    for e in range (0, len(xevents)):
        eventpar = xevents[e].findall('d:origin',ns)
        eventarr = eventpar[0].findall('d:arrival',ns)
        eventpha = eventarr[0].findall('d:phase',ns)
        eventmars = eventpar[0].findall('mars:distance',ns)
        if len(eventmars)>0:
            eventdist = eventmars[0].findall('mars:value', ns)
            eventdescr = xevents[e].findall('d:description',ns)
            check= eventdescr[0].findall('d:type',ns)
            checktxt=check[0].text
            if checktxt == 'earthquake name':
                eventname = eventdescr[0].findall('d:text',ns)
            else:
                eventname = eventdescr[1].findall('d:text',ns)
            eventqual = eventpar[0].findall('mars:locationQuality',ns)
            eventtype = xevents[e].findall('mars:type',ns)

            eventcre = eventpar[0].findall('d:time', ns)
            eventtime = eventcre[0].findall('d:value', ns)

            name = eventname[0].text
            etime = eventtime[0].text
            epdist = eventdist[0].text
            mtype = eventtype[0].text[53:72]
            mqual = eventqual[0].text[63:68]

            line = [name,' ', etime,' ', epdist,' ', mqual,' ', mtype + '\n']
            listline = name + '  ' + etime + '  ' + epdist + '  ' + mqual + '  ' + mtype
            #
            if Avar==1 and mqual=='A':
                if LFvar==1 and mtype=='LOW_FREQUENCY':
                    catbox.insert(e, listline)
                elif BBvar==1 and mtype=='BROADBAND':
                    catbox.insert(e, listline)
                elif HZvar==1 and mtype=='2.4_HZ':
                    catbox.insert(e, listline)
                elif HFvar==1 and mtype=='HIGH_FREQUENCY':
                    catbox.insert(e, listline)
                elif VFvar==1 and mtype=='VERY_HIGH_FREQUENCY':
                    catbox.insert(e, listline)
            elif Bvar==1 and mqual=='B':
                if LFvar==1 and mtype=='LOW_FREQUENCY':
                    catbox.insert(e, listline)
                elif BBvar==1 and mtype=='BROADBAND':
                    catbox.insert(e, listline)
                elif HZvar==1 and mtype=='2.4_HZ':
                    catbox.insert(e, listline)
                elif HFvar==1 and mtype=='HIGH_FREQUENCY':
                    catbox.insert(e, listline)
                elif VFvar==1 and mtype=='VERY_HIGH_FREQUENCY':
                    catbox.insert(e, listline)
            elif Cvar==1 and mqual=='C':
                if LFvar==1 and mtype=='LOW_FREQUENCY':
                    catbox.insert(e, listline)
                elif BBvar==1 and mtype=='BROADBAND':
                    catbox.insert(e, listline)
                elif HZvar==1 and mtype=='2.4_HZ':
                    catbox.insert(e, listline)
                elif HFvar==1 and mtype=='HIGH_FREQUENCY':
                    catbox.insert(e, listline)
                elif VFvar==1 and mtype=='VERY_HIGH_FREQUENCY':
                    catbox.insert(e, listline)
            elif Dvar==1 and mqual=='D':
                if LFvar==1 and mtype=='LOW_FREQUENCY':
                    catbox.insert(e, listline)
                elif BBvar==1 and mtype=='BROADBAND':
                    catbox.insert(e, listline)
                elif HZvar==1 and mtype=='2.4_HZ':
                    catbox.insert(e, listline)
                elif HFvar==1 and mtype=='HIGH_FREQUENCY':
                    catbox.insert(e, listline)
                elif VFvar==1 and mtype=='VERY_HIGH_FREQUENCY':
                    catbox.insert(e, listline)

############################################
# DATA DOWNLOADING AND BASIC PROCESS
############################################

def select_event(self):
    global eventname
    global catbox
    global eventline
    for i in catbox.curselection():
        eventline = catbox.get(i)
        eventname.set(eventline[0:6])

def eventdownload():
    global eventline

    for widget in waveframe2.winfo_children():
        widget.destroy()
    for widget in specframe2.winfo_children():
        widget.destroy()
    minfreq_st.set(0.01)
    maxfreq_st.set(9.99)
    mintime_st.set(-1800)
    maxtime_st.set(5400)
    amplitude_st.set(0)

    import os
    global catbox
    global eventname
    for i in catbox.curselection():
        eventline = catbox.get(i)
    lst = []
    for pos,char in enumerate(eventline):
        if(char == ' '):
            lst.append(pos)
    event = eventline[0:(lst[0])]
    eventtime = eventline[(lst[1]+1):(lst[2])]
    quality = eventline[(lst[5]+1):(lst[6])]
    eventclass = eventline[(lst[7]+1):len(eventline)]
    #
    time_1 = str(eventtime)
    time_2 = time_1
    class_1 = str(eventclass)
    class_2 = class_1
    quality_1 = str(quality)
    q1l = len(quality_1)-3
    quality_2 = quality_1
    from obspy import UTCDateTime
    time = UTCDateTime(time_1)
    starttime = time - 60*30
    endtime = time + 60*90

    fileraw = ('DATA/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/' + str(eventname.get()) + '.mseed')
    filedisp = ('DATA/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/' + str(eventname.get()) + '.mseed')
    filevel = ('DATA/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/' + str(eventname.get()) + '.mseed')
    fileacc = ('DATA/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/' + str(eventname.get()) + '.mseed')

    checkraw = os.path.isfile(fileraw)
    checkdisp = os.path.isfile(filedisp)
    checkvel = os.path.isfile(filevel)
    checkacc = os.path.isfile(fileacc)

    if (checkraw == True and checkdisp == True and checkvel == True and checkacc == True):
        ranval = 1
    else:
        net = 'XB'
        sta = 'ELYSE'
        loc = '02'
        chan = 'BH*'
        from obspy.clients.fdsn import Client
        client = Client('IRIS')
        st = client.get_waveforms(net, sta, loc, chan, starttime, endtime, attach_response = True)
        for tr in st.select(component='U'):
            st.merge(tr)
        timeE1=st[0].stats.starttime;
        timeN1=st[1].stats.starttime;
        timeZ1=st[2].stats.starttime;
        stime=max(timeE1, timeN1, timeZ1)
        timeEe=st[0].stats.endtime;
        timeNe=st[1].stats.endtime;
        timeZe=st[2].stats.endtime;
        etime=min(timeEe, timeNe, timeZe)
        st[0].trim(stime, etime);
        st[1].trim(stime, etime);
        st[2].trim(stime, etime);
        import obspy
        from obspy import read
        from obspy import read_inventory
        inv = obspy.read_inventory('ELYSE.dataless')
        for q in range(1,4):
            if q==1:
                comp=('DISP')
            elif q==2:
                comp=('VEL')
            elif q==3:
                comp=('ACC')
            st_rem1=st.copy()
            pre_filt = [0.005, 0.01, 8, 10] #for 20 Hz data
            st_rem1.remove_response(output = comp, taper_fraction=0.05, pre_filt = pre_filt, inventory = inv);
            sta = inv[0][0]
            azs = []
            dips = []
            trs = []
            channels = ['BHU','BHV','BHW']
            for chn in channels:
                chndata = sta.select(channel=chn)[0]
                azs.append(chndata.azimuth)
                dips.append(chndata.dip)
            from obspy.signal.rotate import rotate2zne
            (z, n, e) = rotate2zne(st_rem1[0], azs[0], dips[0], st_rem1[1], azs[1], dips[1], st_rem1[2], azs[2], dips[2])
            from scipy import signal
            lenz = len(z)
            alp = 5e-2
            window = signal.tukey(len(z), alpha = alp)
            z = z * window
            n = n * window
            e = e * window
            st_new1=st_rem1.copy()
            st_new1[0].data = z;
            st_new1[0].stats.channel = 'BHZ'
            st_new1[1].data = n;
            st_new1[1].stats.channel = 'BHN'
            st_new1[2].data = e;
            st_new1[2].stats.channel = 'BHE'
            path = ('DATA/')
            check = os.path.isdir(path)
            if check == False:
                os.mkdir(path)
            path = ('DATA/' + class_2 + '/')
            check = os.path.isdir(path)
            if check == False:
                os.mkdir(path)
            path = ('DATA/' + class_2 + '/' + quality_2 + '/')
            check = os.path.isdir(path)
            if check == False:
                os.mkdir(path)
            path = ('DATA/' + class_2 + '/' + quality_2 + '/' + event + '/')
            check = os.path.isdir(path)
            if check == False:
                os.mkdir(path)
            filename1 = (path + '/' + event + '_' + comp + '.mseed')
            st_new1.write(filename1, format = 'MSEED')
        #
        targetfile=(path + '/' + event + '.mseed')
        from shutil import copyfile
        #copyfile('fdsnws_msds.mseed', targetfile)
        st.write(targetfile, format='MSEED')

def calc_spectrogram(trace,winlen_sec,overlap):
    import matplotlib.mlab as mlab
    from obspy.signal.util import next_pow_2
    winlen = int(winlen_sec*trace.stats.sampling_rate)
    Fs = trace.stats.sampling_rate
    noverlap = int(winlen*overlap)
    p,f,t = mlab.specgram(trace.data, NFFT=winlen,
                          Fs=Fs, noverlap=noverlap,
                          pad_to=next_pow_2(winlen)*4)
    return p,f,t

def find_nearest(array, value):
    import numpy as np
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def specmaker():
    from obspy import read
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
    from matplotlib.figure import Figure
    import matplotlib.dates as mdates

    global eventname
    global speccomp
    global eventline
    global mintime_st
    global maxtime_st

    lst = []
    for pos,char in enumerate(eventline):
        if(char == ' '):
            lst.append(pos)
    event = eventline[0:(lst[0])]
    eventtime = eventline[(lst[1]+1):(lst[2])]
    quality = eventline[(lst[5]+1):(lst[6])]
    eventclass = eventline[(lst[7]+1):len(eventline)]

    if (speccomp.get() == '1'):
        file = ('DATA/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/' + str(eventname.get()) + '.mseed')
    elif (speccomp.get() == '2'):
        file = ('DATA/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/' + str(eventname.get()) + '_DISP.mseed')
    elif (speccomp.get() == '3'):
        file = ('DATA/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/' + str(eventname.get()) + '_VEL.mseed')
    elif (speccomp.get() == '4'):
        file = ('DATA/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/' + str(eventname.get()) + '_ACC.mseed')

    #set main parameters
    dirty_glitch_remover = True

    winlen_sec = 50.0
    overlap = 0.9
    bp_fmin = 1/5.0
    bp_fmax = 1/3.0

    #set mdates parameters for plotting
    hours = mdates.HourLocator()
    mins = mdates.MinuteLocator(byminute=np.arange(30,60,30))
    cmap='plasma'

    st = read(file)

    if (speccomp.get() == '2' or speccomp.get() == '3' or speccomp.get() == '4'):
        #calculate Z spectrogram------------------------
        tr = st.select(channel='BHZ').copy()[0]
        time = np.linspace(0,tr.stats.npts*tr.stats.delta,tr.stats.npts)
        bhz_p,f,t = calc_spectrogram(tr,winlen_sec,overlap=overlap)
        #N spectrogram------------------------
        tr = st.select(channel='BHN').copy()[0]
        bhn_p,f,t = calc_spectrogram(tr,winlen_sec,overlap=overlap)
        #E spectrogram------------------------
        tr = st.select(channel='BHE').copy()[0]
        bhe_p,f,t = calc_spectrogram(tr,winlen_sec,overlap=overlap)
    elif (speccomp.get() == '1'):
        #calculate Z spectrogram------------------------
        tr = st.select(channel='BHU').copy()[0]
        time = np.linspace(0,tr.stats.npts*tr.stats.delta,tr.stats.npts)
        bhz_p,f,t = calc_spectrogram(tr,winlen_sec,overlap=overlap)
        #N spectrogram------------------------
        tr = st.select(channel='BHV').copy()[0]
        bhn_p,f,t = calc_spectrogram(tr,winlen_sec,overlap=overlap)
        #E spectrogram------------------------
        tr = st.select(channel='BHW').copy()[0]
        bhe_p,f,t = calc_spectrogram(tr,winlen_sec,overlap=overlap)

    for widget in specframe2.winfo_children():
        widget.destroy()

    fig = plt.figure(dpi=70, layout='constrained')
    plt.clf()
    cm=2.54
    fig.set_size_inches(164/cm, 164/cm)
    gs = fig.add_gridspec(1, 3)
    axs = gs.subplots(sharex = False, sharey=False)

    if (speccomp.get() == '1'):
        vmin = -10
        vmax = 50
    if (speccomp.get() == '2'):
        vmin = -250
        vmax = -200
    elif (speccomp.get() == '3'):
        vmin = -230 #dB (for colorbar)
        vmax = -180 #dB (for colorbar)
    elif (speccomp.get() == '4'):
        vmin = -210
        vmax = -160

    tw = [-1800, 5400];
    lfw = [minfreq_st.get(), minfreq_st.get()]
    hfw = [maxfreq_st.get(), maxfreq_st.get()]

    fw = [0, 10];
    ltw = [mintime_st.get(), mintime_st.get()]
    htw = [maxtime_st.get(), maxtime_st.get()]

    t = t-1800;

    ###
    axs[0].pcolormesh(t,f,10*np.log10(bhz_p), vmin=vmin,vmax=vmax,cmap=cmap);
    axs[0].plot(tw,lfw, '--', color = 'white');
    axs[0].plot(tw,hfw, '--', color = 'white');
    axs[0].plot(ltw,fw, '--', color = 'white');
    axs[0].plot(htw,fw, '--', color = 'white');
    titlestrZ='BHZ'
    if (speccomp.get() == '1'):
        titlestrZ = 'BHU'
    axs[0].set_title(titlestrZ);
    axs[0].set_ylabel('Frequency (Hz)');
    axs[0].set_xlabel('Time (s)');
    axs[0].set_xlim(-1800, 5400)
    axs[0].patch.set_facecolor('#ececec')
    #
    axs[1].pcolormesh(t,f,10*np.log10(bhn_p), vmin=vmin,vmax=vmax,cmap=cmap);
    axs[1].plot(tw,lfw, '--', color = 'white');
    axs[1].plot(tw,hfw, '--', color = 'white');
    axs[1].plot(ltw,fw, '--', color = 'white');
    axs[1].plot(htw,fw, '--', color = 'white');
    titlestrN='BHN'
    if (speccomp.get() == '1'):
        titlestrN = 'BHV'
    axs[1].set_title(titlestrN);
    axs[1].set_ylabel('Frequency (Hz)');
    axs[1].set_xlabel('Time (s)');
    axs[1].set_xlim(-1800, 5400)
    axs[1].patch.set_facecolor('#ececec')
    #
    axs[2].pcolormesh(t,f,10*np.log10(bhn_p), vmin=vmin,vmax=vmax,cmap=cmap);
    axs[2].plot(tw,lfw, '--', color = 'white');
    axs[2].plot(tw,hfw, '--', color = 'white');
    axs[2].plot(ltw,fw, '--', color = 'white');
    axs[2].plot(htw,fw, '--', color = 'white');
    titlestrE='BHE'
    if (speccomp.get() == '1'):
        titlestrE = 'BHW'
    axs[2].set_title(titlestrE);
    axs[2].set_ylabel('Frequency (Hz)');
    axs[2].set_xlabel('Time (s)');
    axs[2].set_xlim(-1800, 5400)
    axs[2].patch.set_facecolor('#ececec')


    if (speccomp.get() == '1'):
        suptitletext = str(eventname.get()) + ' | Raw data | Event time: ' + eventtime;
    elif (speccomp.get() == '2'):
        suptitletext = str(eventname.get()) + ' | Displacement | Event time: ' + eventtime;
    elif (speccomp.get() == '3'):
        suptitletext = str(eventname.get()) + ' | Velocity | Event time: ' + eventtime;
    elif (speccomp.get() == '4'):
        suptitletext = str(eventname.get()) + ' | Acceleration | Event time: ' + eventtime;

    fig.suptitle(suptitletext)

    fig.patch.set_facecolor('#ececec')
    #fig2.tight_layout()
    canvaswaves = FigureCanvasTkAgg(fig, specframe2)
    canvaswaves.draw()
    canvaswaves.get_tk_widget().pack()

def specfigmaker():

    from obspy import read
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
    from matplotlib.figure import Figure
    import matplotlib.dates as mdates

    global eventname
    global speccomp
    global eventline
    global mintime_st
    global maxtime_st

    lst = []
    for pos,char in enumerate(eventline):
        if(char == ' '):
            lst.append(pos)
    event = eventline[0:(lst[0])]
    eventtime = eventline[(lst[1]+1):(lst[2])]
    quality = eventline[(lst[5]+1):(lst[6])]
    eventclass = eventline[(lst[7]+1):len(eventline)]

    if (speccomp.get() == '1'):
        file = ('DATA/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/' + str(eventname.get()) + '.mseed')
    elif (speccomp.get() == '2'):
        file = ('DATA/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/' + str(eventname.get()) + '_DISP.mseed')
    elif (speccomp.get() == '3'):
        file = ('DATA/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/' + str(eventname.get()) + '_VEL.mseed')
    elif (speccomp.get() == '4'):
        file = ('DATA/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/' + str(eventname.get()) + '_ACC.mseed')

    #set main parameters
    dirty_glitch_remover = True

    winlen_sec = 50.0
    overlap = 0.9
    bp_fmin = 1/5.0
    bp_fmax = 1/3.0

    #set mdates parameters for plotting
    hours = mdates.HourLocator()
    mins = mdates.MinuteLocator(byminute=np.arange(30,60,30))
    cmap='plasma'

    st = read(file)

    if (speccomp.get() == '2' or speccomp.get() == '3' or speccomp.get() == '4'):
        #calculate Z spectrogram------------------------
        tr = st.select(channel='BHZ').copy()[0]
        time = np.linspace(0,tr.stats.npts*tr.stats.delta,tr.stats.npts)
        bhz_p,f,t = calc_spectrogram(tr,winlen_sec,overlap=overlap)
        #N spectrogram------------------------
        tr = st.select(channel='BHN').copy()[0]
        bhn_p,f,t = calc_spectrogram(tr,winlen_sec,overlap=overlap)
        #E spectrogram------------------------
        tr = st.select(channel='BHE').copy()[0]
        bhe_p,f,t = calc_spectrogram(tr,winlen_sec,overlap=overlap)
    elif (speccomp.get() == '1'):
        #calculate Z spectrogram------------------------
        tr = st.select(channel='BHU').copy()[0]
        time = np.linspace(0,tr.stats.npts*tr.stats.delta,tr.stats.npts)
        bhz_p,f,t = calc_spectrogram(tr,winlen_sec,overlap=overlap)
        #N spectrogram------------------------
        tr = st.select(channel='BHV').copy()[0]
        bhn_p,f,t = calc_spectrogram(tr,winlen_sec,overlap=overlap)
        #E spectrogram------------------------
        tr = st.select(channel='BHW').copy()[0]
        bhe_p,f,t = calc_spectrogram(tr,winlen_sec,overlap=overlap)

    fig3 = plt.figure(dpi=200, layout='constrained')
    fig3.subplots_adjust(wspace=0.5)
    plt.clf()
    cm=2.54
    fig3.set_size_inches(40/cm, 10/cm)
    gs = fig3.add_gridspec(1, 3)
    axs = gs.subplots(sharex = False, sharey=False)

    if (speccomp.get() == '1'):
        vmin = -10
        vmax = 50
    if (speccomp.get() == '2'):
        vmin = -250
        vmax = -200
    elif (speccomp.get() == '3'):
        vmin = -230 #dB (for colorbar)
        vmax = -180 #dB (for colorbar)
    elif (speccomp.get() == '4'):
        vmin = -210
        vmax = -160

    tw = [-1800, 5400];
    lfw = [minfreq_st.get(), minfreq_st.get()]
    hfw = [maxfreq_st.get(), maxfreq_st.get()]

    fw = [0, 10];
    ltw = [mintime_st.get(), mintime_st.get()]
    htw = [maxtime_st.get(), maxtime_st.get()]

    t = t-1800;

    ###
    axs[0].pcolormesh(t,f,10*np.log10(bhz_p), vmin=vmin,vmax=vmax,cmap=cmap);
    axs[0].plot(tw,lfw, '--', color = 'white');
    axs[0].plot(tw,hfw, '--', color = 'white');
    axs[0].plot(ltw,fw, '--', color = 'white');
    axs[0].plot(htw,fw, '--', color = 'white');
    titlestrZ='BHZ'
    if (speccomp.get() == '1'):
        titlestrZ = 'BHU'
    axs[0].set_title(titlestrZ);
    axs[0].set_ylabel('Frequency (Hz)');
    axs[0].set_xlabel('Time (s)');
    axs[0].set_xlim(-1800, 5400)
    axs[0].patch.set_facecolor('#ffffff')
    #
    axs[1].pcolormesh(t,f,10*np.log10(bhn_p), vmin=vmin,vmax=vmax,cmap=cmap);
    axs[1].plot(tw,lfw, '--', color = 'white');
    axs[1].plot(tw,hfw, '--', color = 'white');
    axs[1].plot(ltw,fw, '--', color = 'white');
    axs[1].plot(htw,fw, '--', color = 'white');
    titlestrN='BHN'
    if (speccomp.get() == '1'):
        titlestrN = 'BHV'
    axs[1].set_title(titlestrN);
    axs[1].set_ylabel('Frequency (Hz)');
    axs[1].set_xlabel('Time (s)');
    axs[1].set_xlim(-1800, 5400)
    axs[1].patch.set_facecolor('#ffffff')
    #
    axs[2].pcolormesh(t,f,10*np.log10(bhn_p), vmin=vmin,vmax=vmax,cmap=cmap);
    axs[2].plot(tw,lfw, '--', color = 'white');
    axs[2].plot(tw,hfw, '--', color = 'white');
    axs[2].plot(ltw,fw, '--', color = 'white');
    axs[2].plot(htw,fw, '--', color = 'white');
    titlestrE='BHE'
    if (speccomp.get() == '1'):
        titlestrE = 'BHW'
    axs[2].set_title(titlestrE);
    axs[2].set_ylabel('Frequency (Hz)');
    axs[2].set_xlabel('Time (s)');
    axs[2].set_xlim(-1800, 5400)
    axs[2].patch.set_facecolor('#ffffff')

    if (speccomp.get() == '1'):
        suptitletext = str(eventname.get()) + ' | Raw data | Event time: ' + eventtime;
    elif (speccomp.get() == '2'):
        suptitletext = str(eventname.get()) + ' | Displacement | Event time: ' + eventtime;
    elif (speccomp.get() == '3'):
        suptitletext = str(eventname.get()) + ' | Velocity | Event time: ' + eventtime;
    elif (speccomp.get() == '4'):
        suptitletext = str(eventname.get()) + ' | Acceleration | Event time: ' + eventtime;

    fig3.suptitle(suptitletext)

    fig3.patch.set_facecolor('#ffffff')

    import os
    path = ('Figures/')
    check = os.path.isdir(path)
    if check == False:
        os.mkdir(path)
    path = ('Figures/Spectrograms/')
    check = os.path.isdir(path)
    if check == False:
        os.mkdir(path)
    path = ('Figures/Spectrograms/' + str(eventclass) + '/')
    check = os.path.isdir(path)
    if check == False:
        os.mkdir(path)
    path = ('Figures/Spectrograms/' + str(eventclass) + '/' + str(quality) + '/')
    check = os.path.isdir(path)
    if check == False:
        os.mkdir(path)
    path = ('Figures/Spectrograms/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/')
    check = os.path.isdir(path)
    if check == False:
        os.mkdir(path)

    if (speccomp.get() == '1'):
        specfigfile = (path + '/' + str(eventname.get()) + '_raw.png')
    elif  (speccomp.get() == '2'):
        specfigfile = (path + '/' + str(eventname.get()) + '_DISP.png')
    elif (speccomp.get() == '3'):
        specfigfile = (path + '/' + str(eventname.get()) + '_VEL.png')
    elif (speccomp.get() == '4'):
        specfigfile = (path + '/' + str(eventname.get()) + '_ACC.png')

    fig3.savefig(specfigfile)

def wavemaker():
    from obspy import read
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
    from matplotlib.figure import Figure

    global eventname
    global speccomp
    global eventline
    global minfreq_st
    global maxfreq_st
    global mintime_st
    global maxtime_st
    global amplitude_st
    global eventline

    for widget in waveframe2.winfo_children():
        widget.destroy()

    lst = []
    for pos,char in enumerate(eventline):
        if(char == ' '):
            lst.append(pos)
    event = eventline[0:(lst[0])]
    eventtime = eventline[(lst[1]+1):(lst[2])]
    quality = eventline[(lst[5]+1):(lst[6])]
    eventclass = eventline[(lst[7]+1):len(eventline)]

    specmaker()

    fileraw = ('DATA/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/' + str(eventname.get()) + '.mseed')
    filedisp = ('DATA/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/' + str(eventname.get()) + '_DISP.mseed')
    filevel = ('DATA/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/' + str(eventname.get()) + '_VEL.mseed')
    fileacc = ('DATA/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/' + str(eventname.get()) + '_ACC.mseed')

    if (speccomp.get() == '1'):
        st = read(fileraw)
    elif  (speccomp.get() == '2'):
        st = read(filedisp)
    elif (speccomp.get() == '3'):
        st = read(filevel)
    elif (speccomp.get() == '4'):
        st = read(fileacc)

    timevect = st[0].times() - 1800;
    Zdata = st[0].data
    Ndata = st[1].data
    Edata = st[2].data

    from scipy import signal
    sf = timevect[1]-timevect[0];
    fs=1/sf;
    nyq = 0.5 * fs;
    lowcut = float(minfreq_st.get());
    highcut = float(maxfreq_st.get());
    low = lowcut / nyq;
    high = highcut / nyq;
    b, a = signal.butter(2, [low, high], btype='band')
    Zf = signal.filtfilt(b, a, Zdata)
    Nf = signal.filtfilt(b, a, Ndata)
    Ef = signal.filtfilt(b, a, Edata)

    mint = float(mintime_st.get());
    maxt = float(maxtime_st.get());

    fig2 = plt.figure(dpi=70, layout='constrained')
    plt.clf()
    cm=2.54
    fig2.set_size_inches(164/cm, 164/cm)
    gs = fig2.add_gridspec(3, 1)
    axs = gs.subplots(sharex = True, sharey=True)

    if (speccomp.get() == '1'):
        yvalues_txt = 'DU'
    elif  (speccomp.get() == '2'):
        yvalues_txt = 'Displacement (m)'
    elif (speccomp.get() == '3'):
        yvalues_txt = 'Velocity (m/s)'
    elif (speccomp.get() == '4'):
        yvalues_txt = 'Acceleration (m/s^2)'

    # Now trim

    start_t = find_nearest(timevect, mint)
    stop_t = find_nearest(timevect, maxt)

    timef = timevect[start_t:stop_t]
    Zfn = Zf[start_t:stop_t]
    Nfn = Nf[start_t:stop_t]
    Efn = Ef[start_t:stop_t]


    ###
    axs[0].plot(timef,Zfn,'k-', linewidth=0.5);
    titlestrZ='BHZ | f = ' + str(lowcut) + '-' + str(highcut) + ' Hz'
    if (speccomp.get() == '1'):
        titlestrZ='BHU | f = ' + str(lowcut) + '-' + str(highcut) + ' Hz'
    axs[0].set_title(titlestrZ);
    axs[0].set_ylabel(yvalues_txt);
    axs[0].set_xlabel('Time (s)');
    axs[0].set_xlim(mint, maxt);
    if float(amplitude_st.get()) > 0:
        axs[0].set_ylim(-float(amplitude_st.get()), float(amplitude_st.get()));
    axs[0].patch.set_facecolor('#ececec')
    #
    axs[1].plot(timef,Nfn,'k-', linewidth=0.5);
    titlestrN='BHN | f = ' + str(lowcut) + '-' + str(highcut) + ' Hz'
    if (speccomp.get() == '1'):
        titlestrN='BHV | f = ' + str(lowcut) + '-' + str(highcut) + ' Hz'
    axs[1].set_title(titlestrN);
    axs[1].set_ylabel(yvalues_txt);
    axs[1].set_xlabel('Time (s)');
    axs[1].set_xlim(mint, maxt);
    axs[1].patch.set_facecolor('#ececec')
    #
    axs[2].plot(timef,Efn,'k-', linewidth=0.5);
    titlestrE='BHE | f = ' + str(lowcut) + '-' + str(highcut) + ' Hz'
    if (speccomp.get() == '1'):
        titlestrE='BHW | f = ' + str(lowcut) + '-' + str(highcut) + ' Hz'
    axs[2].set_title(titlestrE);
    axs[2].set_ylabel(yvalues_txt);
    axs[2].set_xlabel('Time (s)');
    axs[2].set_xlim(mint, maxt);
    axs[2].patch.set_facecolor('#ececec')

    suptitletext = str(eventname.get()) + ' | Event time: ' + eventtime;
    fig2.suptitle(suptitletext)

    fig2.patch.set_facecolor('#ececec')
    #fig2.tight_layout()
    canvaswaves = FigureCanvasTkAgg(fig2, waveframe2)
    canvaswaves.draw()
    canvaswaves.get_tk_widget().pack()

    print ('It works fine')

def seisfigmaker():

    from obspy import read
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
    from matplotlib.figure import Figure

    global eventname
    global speccomp
    global eventline
    global minfreq_st
    global maxfreq_st
    global mintime_st
    global maxtime_st
    global amplitude_st
    global eventline

    lst = []
    for pos,char in enumerate(eventline):
        if(char == ' '):
            lst.append(pos)
    event = eventline[0:(lst[0])]
    eventtime = eventline[(lst[1]+1):(lst[2])]
    quality = eventline[(lst[5]+1):(lst[6])]
    eventclass = eventline[(lst[7]+1):len(eventline)]

    fileraw = ('DATA/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/' + str(eventname.get()) + '.mseed')
    filedisp = ('DATA/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/' + str(eventname.get()) + '_DISP.mseed')
    filevel = ('DATA/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/' + str(eventname.get()) + '_VEL.mseed')
    fileacc = ('DATA/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/' + str(eventname.get()) + '_ACC.mseed')

    if (speccomp.get() == '1'):
        st = read(fileraw)
    elif  (speccomp.get() == '2'):
        st = read(filedisp)
    elif (speccomp.get() == '3'):
        st = read(filevel)
    elif (speccomp.get() == '4'):
        st = read(fileacc)

    timevect = st[0].times() - 1800;
    Zdata = st[0].data
    Ndata = st[1].data
    Edata = st[2].data

    from scipy import signal
    sf = timevect[1]-timevect[0];
    fs=1/sf;
    nyq = 0.5 * fs;
    lowcut = float(minfreq_st.get());
    highcut = float(maxfreq_st.get());
    low = lowcut / nyq;
    high = highcut / nyq;
    b, a = signal.butter(2, [low, high], btype='band')
    Zf = signal.filtfilt(b, a, Zdata)
    Nf = signal.filtfilt(b, a, Ndata)
    Ef = signal.filtfilt(b, a, Edata)

    mint = float(mintime_st.get());
    maxt = float(maxtime_st.get());

    fig4 = plt.figure(dpi=200, layout='constrained')
    fig4.subplots_adjust(hspace=0.5)

    plt.clf()
    cm=2.54
    fig4.set_size_inches(40/cm, 20/cm)
    gs = fig4.add_gridspec(3, 1)
    axs = gs.subplots(sharex = True, sharey=True)

    if (speccomp.get() == '1'):
        yvalues_txt = 'DU'
    elif  (speccomp.get() == '2'):
        yvalues_txt = 'Displacement (m)'
    elif (speccomp.get() == '3'):
        yvalues_txt = 'Velocity (m/s)'
    elif (speccomp.get() == '4'):
        yvalues_txt = 'Acceleration (m/s^2)'

    # Now trim

    start_t = find_nearest(timevect, mint)
    stop_t = find_nearest(timevect, maxt)

    timef = timevect[start_t:stop_t]
    Zfn = Zf[start_t:stop_t]
    Nfn = Nf[start_t:stop_t]
    Efn = Ef[start_t:stop_t]


    ###
    axs[0].plot(timef,Zfn,'k-', linewidth=0.5);
    titlestrZ='BHZ | f = ' + str(lowcut) + '-' + str(highcut) + ' Hz'
    if (speccomp.get() == '1'):
        titlestrZ='BHU | f = ' + str(lowcut) + '-' + str(highcut) + ' Hz'
    axs[0].set_title(titlestrZ);
    axs[0].set_ylabel(yvalues_txt);
    axs[0].set_xlabel('Time (s)');
    axs[0].set_xlim(mint, maxt);
    if float(amplitude_st.get()) > 0:
        axs[0].set_ylim(-float(amplitude_st.get()), float(amplitude_st.get()));
    axs[0].patch.set_facecolor('#ffffff')
    #
    axs[1].plot(timef,Nfn,'k-', linewidth=0.5);
    titlestrN='BHN | f = ' + str(lowcut) + '-' + str(highcut) + ' Hz'
    if (speccomp.get() == '1'):
        titlestrN='BHV | f = ' + str(lowcut) + '-' + str(highcut) + ' Hz'
    axs[1].set_title(titlestrN);
    axs[1].set_ylabel(yvalues_txt);
    axs[1].set_xlabel('Time (s)');
    axs[1].set_xlim(mint, maxt);
    axs[1].patch.set_facecolor('#ffffff')
    #
    axs[2].plot(timef,Efn,'k-', linewidth=0.5);
    titlestrE='BHE | f = ' + str(lowcut) + '-' + str(highcut) + ' Hz'
    if (speccomp.get() == '1'):
        titlestrE='BHW | f = ' + str(lowcut) + '-' + str(highcut) + ' Hz'
    axs[2].set_title(titlestrE);
    axs[2].set_ylabel(yvalues_txt);
    axs[2].set_xlabel('Time (s)');
    axs[2].set_xlim(mint, maxt);
    axs[2].patch.set_facecolor('#ffffff')

    fig4.patch.set_facecolor('#ffffff')

    suptitletext = str(eventname.get()) + ' | Event time: ' + eventtime;
    fig4.suptitle(suptitletext)

    print ('It works fine')

    import os
    path = ('Figures/')
    check = os.path.isdir(path)
    if check == False:
        os.mkdir(path)
    path = ('Figures/Seismograms/')
    check = os.path.isdir(path)
    if check == False:
        os.mkdir(path)
    path = ('Figures/Seismograms/' + str(eventclass) + '/')
    check = os.path.isdir(path)
    if check == False:
        os.mkdir(path)
    path = ('Figures/Seismograms/' + str(eventclass) + '/' + str(quality) + '/')
    check = os.path.isdir(path)
    if check == False:
        os.mkdir(path)
    path = ('Figures/Seismograms/' + str(eventclass) + '/' + str(quality) + '/' + str(eventname.get()) + '/')
    check = os.path.isdir(path)
    if check == False:
        os.mkdir(path)

    if (speccomp.get() == '1'):
        seisfigfile = (path + '/' + str(eventname.get()) + '_raw.png')
    elif  (speccomp.get() == '2'):
        seisfigfile = (path + '/' + str(eventname.get()) + '_DISP.png')
    elif (speccomp.get() == '3'):
        seisfigfile = (path + '/' + str(eventname.get()) + '_VEL.png')
    elif (speccomp.get() == '4'):
        seisfigfile = (path + '/' + str(eventname.get()) + '_ACC.png')

    fig4.savefig(seisfigfile)

def our_command():
    pass

################################################################################
# Make the interface
################################################################################

# Create menu items

my_menu = Menu(root)

root.config(menu = my_menu)

file_menu = Menu(my_menu)
my_menu.add_cascade(label='File', menu=file_menu)
file_menu.add_command(label='New...', command=our_command)

help_menu = Menu(my_menu)
my_menu.add_cascade(label='Help', menu = help_menu)
help_menu.add_command(label='New...', command=our_command)

# Add the catalog xml file

# Put the catalog on the left side of the Interface
catframe = Frame(root, width = 345, height = 740, borderwidth=1, relief='solid')
catframe.pack_propagate(0)
catframe.place(x=5, y=5)

# Button for catalog creation
runButton = Button(catframe, text='SEARCH EVENTS IN CATALOG', fg='black', font='Arial 14 bold', command=makecatalog)
runButton.config(width=37, height=2)
runButton.place(x=5, y=5)
runButton.pack_propagate(0)

catframe2 = Frame(catframe, width = 335, height = 680, borderwidth=0, relief='solid')
catframe2.pack_propagate(0)
catframe2.place(x=5,y=50)
# Add the options
qualitytxt = Label(catframe2, text='Quality', font = 'Arial 11 bold')
qualitytxt.place(x=12, y=0)

Amode = IntVar()
Abutton = Checkbutton(catframe2, text='A', variable=Amode, onvalue = 1, offvalue = 0, font = 'Arial 10', state='active')
Abutton.place(x=0, y=20)
#
Bmode = IntVar()
Bbutton = Checkbutton(catframe2, text='B', variable=Bmode, onvalue = 1, offvalue = 0, font = 'Arial 10', state='active')
Bbutton.place(x=35, y=20)
#
Cmode = IntVar()
Cbutton = Checkbutton(catframe2, text='C', variable=Cmode, onvalue = 1, offvalue = 0, font = 'Arial 10', state='active')
Cbutton.place(x=0, y=40)
#
Dmode = IntVar()
Dbutton = Checkbutton(catframe2, text='D', variable=Dmode, onvalue = 1, offvalue = 0, font = 'Arial 10', state='active')
Dbutton.place(x=35, y=40)

canvas=Canvas(catframe2, width=10, height=50)
canvas.place(x=70, y=10)
canvas.create_line(5,5,5,50, fill='black', width=1)

freqtypetxt = Label(catframe2, text='Frequency Type', font = 'Arial 11 bold')
freqtypetxt.place(x=160, y=0)

LFmode = IntVar()
LFbutton = Checkbutton(catframe2, text = 'Low Frequency', variable = LFmode, onvalue=1, offvalue = 0, font = 'Arial 10', state = 'active')
LFbutton.place(x=80, y=20)
#
BBmode = IntVar()
BBbutton = Checkbutton(catframe2, text = 'Broadband', variable = BBmode, onvalue = 1, offvalue = 0, font = 'Arial 10', state = 'active')
BBbutton.place(x=180, y=20)
#
HZmode = IntVar()
HZbutton = Checkbutton(catframe2, text = '2.4 Hz', variable = HZmode, onvalue = 1, offvalue = 0, font = 'Arial 10', state = 'active')
HZbutton.place(x=265, y=20)
#
HFmode = IntVar()
HFbutton = Checkbutton(catframe2, text = 'High Frequency', variable = HFmode, onvalue = 1, offvalue = 0, font = 'Arial 10', state = 'active')
HFbutton.place(x=80, y=40)
#
VFmode = IntVar()
VFbutton = Checkbutton(catframe2, text = 'Very High Frequency', variable = VFmode, onvalue = 1, offvalue = 0, font = 'Arial 10', state = 'active')
VFbutton.place(x=180, y=40)

eventname = StringVar()
# Add the list box for the catalog
events = StringVar()
catbox = Listbox(catframe2, font = 'Arial 10', borderwidth=0, listvariable = events, selectmode=SINGLE)
catbox.bind('<<ListboxSelect>>', select_event)
catbox.config(width=56, height = 48)
catbox.place(x = -1, y=62)
catscroll = Scrollbar(catframe2, orient='vertical')
#catscroll.config(width = 20)
#catscroll.place(x = 325, y=5)
catbox.config(yscrollcommand = catscroll.set)
catscroll.config(command = catbox.yview)

eventtxt = Label(catframe2, width = 20, text = 'Selected event', font = 'Arial 12')
eventtxt.place(x=5, y=638)
#
event_e = Entry(catframe2, width = 20, font = 'Arial 12', textvariable=eventname, justify = 'center')
event_e.place(x=5, y=656)
#
eventbutton = Button(catframe2, width=20, height=2, text = 'Download and process data', font = 'Arial 12', command=eventdownload)
eventbutton.place(x=158, y=643)

##########################

specframe = Frame(root, width = 840, height = 340, borderwidth=1, relief='solid')
specframe.place(x=355, y=5)

#
speccomp = StringVar(value = 1)
compspecbutton1 = Radiobutton(specframe, text = 'Instrument Counts', variable = speccomp, value = 1, font = 'Arial 11')
compspecbutton2 = Radiobutton(specframe, text = 'Displacement', variable = speccomp, value = 2, font = 'Arial 11')
compspecbutton3 = Radiobutton(specframe, text = 'Velocity', variable = speccomp, value = 3, font = 'Arial 11')
compspecbutton4 = Radiobutton(specframe, text = 'Acceleration', variable = speccomp, value = 4, font = 'Arial 11')

compspecbutton1.place(x=5, y=5)
compspecbutton2.place(x=130, y=5)
compspecbutton3.place(x=235, y=5)
compspecbutton4.place(x=310, y=5)

specbutton = Button(specframe, text = 'Make the spectrograms', width = 14, height = 1, font = 'Arial 11', command=specmaker)
specbutton.place(x=595, y=3)

specfigbutton = Button(specframe, text = 'Export the figure', width = 10, height = 1, font = 'Arial 11', command=specfigmaker)
specfigbutton.place(x=730, y=3)

specframe2 = Frame(specframe, width = 825, height = 285, borderwidth=0, relief = 'solid')
specframe2.place(x=5, y=45)
specframe2.pack_propagate(0)

#

waveframe = Frame(root, width = 840, height =395, borderwidth=1, relief='solid')
waveframe.place(x=355, y=350)

minfreq_st = DoubleVar()
minfreq_st.set(0.01)
maxfreq_st = DoubleVar()
maxfreq_st.set(9.99)

mintime_st = DoubleVar()
mintime_st.set(-1800)
maxtime_st = DoubleVar()
maxtime_st.set(5400)

amplitude_st = DoubleVar()
amplitude_st.set(0)

minfreqtxt = Label(waveframe, text = 'Min freq (Hz)', font = 'Arial 11')
minfreqtxt.place(x=5, y=5)
maxfreqtxt = Label(waveframe, text = 'Max freq (Hz)', font = 'Arial 11')
maxfreqtxt.place(x=115, y=5)

minfreq_e = Entry(waveframe, width = 4, font = 'Arial 11', textvariable=minfreq_st)
minfreq_e.place(x=75, y=3)
maxfreq_e = Entry(waveframe, width = 4, font = 'Arial 11', textvariable=maxfreq_st)
maxfreq_e.place(x=185, y=3)

mintimetxt = Label(waveframe, text = 'Min time (s)', font = 'Arial 11')
mintimetxt.place(x=225, y=5)
maxtimetxt = Label(waveframe, text = 'Max time (s)', font = 'Arial 11')
maxtimetxt.place(x=340, y=5)

mintime_e = Entry(waveframe, width = 5, font = 'Arial 11', textvariable=mintime_st)
mintime_e.place(x=290, y=3)
maxtime_e = Entry(waveframe, width = 5, font = 'Arial 11', textvariable=maxtime_st)
maxtime_e.place(x=405, y=3)

amplitudetxt = Label(waveframe, text = 'Amplitude', font = 'Arial 11')
amplitudetxt.place(x=455, y=5)

amplitude_e = Entry(waveframe, width = 6, font = 'Artial 11', textvariable=amplitude_st)
amplitude_e.place(x=510, y=3)

wavebutton = Button(waveframe, text = 'Show the seismograms', width = 14, height = 1, font = 'Arial 11', command=wavemaker)
wavebutton.place(x=595, y=2)

figseisbutton = Button(waveframe, text = 'Export the figure', width = 10, height = 1, font = 'Arial 11', command=seisfigmaker)
figseisbutton.place(x=730, y=2)

waveframe2 = Frame(waveframe, width = 830, height = 360, borderwidth=0, relief='solid')
waveframe2.place(x=5, y=30)
waveframe2.pack_propagate(0)


root.mainloop()
