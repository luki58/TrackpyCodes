 #!/usr/bin/env python
# coding: utf-8


# Dieses Skript dient zum einen dazu, die richtige Parameterwahl für das Finden und Verbinden von Partikeln zu finden und zum 
# anderen dazu, zu entscheiden, ob der Partikeldrift per trackpy oder manuell berechnet werden soll. Eine manuelle Berechnung 
# ist notwendig, wenn es zu wenige driftende Partikel gibt (trackpy berechnet dann einen zu großen Wert für den Drift).

from __future__ import division, unicode_literals, print_function # Für die Kompatibilität mit Python 2 und 3.
import time
startzeit = time.time()
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pims
import trackpy as tp
import os
from PIL import Image
import json

# Ändere das folgende zu %matplotlib notebook für interaktives plotting.
get_ipython().run_line_magic('matplotlib', 'inline')

# Style anpassen:
mpl.rc('figure', figsize=(19.2, 10.8))
mpl.rc('image', cmap='gray')



path_prefix = 'D:/Backup/Uni/Promotion/Forschung/PK-4/Laser_induced_Shear_Flow/Laser_Data/';
json_name = 'laser_parameter.json';


with open(path_prefix + json_name) as f:
  j_data = json.load(f)


laser_scan = j_data['scanname'];
path_folder = j_data[laser_scan]['Foldername blabla'];


xiq_path = path_prefix + path_folder;
crop_path = path_prefix + path_folder + 'cropped/';
results_path = path_prefix + path_folder + 'Ergebnisse/';

crop_mode = j_data[laser_scan]['cropmode'];

init_frame = j_data[laser_scan]['initframe'];
frame_anz = j_data[laser_scan]['framenumber'];
final_frame = init_frame + frame_anz;

start=init_frame;
stop=final_frame;

# Festlegen des Blidausschnitts
# grob abgeschätzte Mitte des Laserstrahls: 894 Px
# 

y_laser = j_data[laser_scan]['laserpos'];

cropsize_x = j_data[laser_scan]['cropsize_x'];
cropsize_y = j_data[laser_scan]['cropsize_y'];


crop_x2 = j_data[laser_scan]['crop_x2'];
crop_x1 = crop_x2-cropsize_x;

crop_y1 = y_laser-(cropsize_y/2);
crop_y2 = y_laser+(cropsize_y/2);


cut_step = j_data[laser_scan]['cutstep'];

frame_rate = j_data[laser_scan]['framerate'];
spatial_resolution_mu_to_pix = j_data[laser_scan]['mum_2_pixel'];    # mu/pix
spatial_resolution_mm_to_pix = j_data[laser_scan]['millim_2_pixel'];    # mm/pix
scale_v_mu_to_pix = spatial_resolution_mu_to_pix*frame_rate; # mu/pix/s

if not os.path.exists(crop_path):
	os.makedirs(crop_path);

if not os.path.exists(results_path):
	os.makedirs(results_path);


# Funktion zum Zuschneiden der Bilder.
# Wir betrachten nur die rechte Bildhälfte, da sie besser beleuchtet ist.


def crop(x_path, c_path, step, x1, x2, y1, y2, init_f, final_f):
	listdir = sorted(os.listdir(x_path));
	bl = len(str(len(listdir)));
	index = init_f;
	cnt = 0;

	while index < final_f:
		image_obj = Image.open(x_path+'/'+listdir[index]);
		cropped_image = image_obj.crop( (x1, y1, x2, y2) );
		pic_name = 'crop'+(bl-len(str(cnt+1)))*'0'+str(cnt+1)+'.png';
		cropped_image.save(c_path+'/'+pic_name);

		index = index + step;
		cnt = cnt + 1;

	return 0;


if crop_mode==1:
    print('Cropping pictures');
    # Ausschneiden der Bildausschnitte
    listdir = os.listdir(crop_path);

    if len(listdir)!=0:
        for i in range(len(listdir)):
            os.remove(crop_path+listdir[i]);
    # Ausschneiden der Bildausschnitte
    crop(xiq_path, crop_path, cut_step, crop_x1, crop_x2, crop_y1, crop_y2, init_frame,final_frame);


# Festlegen der Bildsequenz
frames = pims.ImageSequence(crop_path +'/*.png', as_grey=True);




#parabola_number = '7' # Parabel-Nr.
#time = '16s' # Das betrachtete Zeitintervall (der angegebene Zeitpunkt ist der exakte Zeitpunkt des ersten Bildes).
# Umwandlung der Farbbilder in Graustufenbilder und Zuschneiden der Bilder. Mit frames werden die Bilder aus dem angegebenen 
# Ordner aufgerufen.


 # Anzahl der betrachteten Bilder.
frames # Zur Überprüfung




print(frames[0].shape) # Zur Überprüfung. Die Umwandlung in Graustufen und das Zuschneiden finden erst statt, wenn ein 
                       # bestimmtes Bild über frames[#] aufgerufen wird.



plt.figure(figsize=(14,14))
plt.imshow(frames[0]); # Erstes frame mit Koordinatenachsen anzeigen. Zur Überprüfung. 



# Bestimmung der Positionen der Mikropartikel auf frames[first_frame] 
L = j_data[laser_scan]['featurelength'] # Länge der größten Linie. Zahl muss ungerade sein! Großzügig abschätzen.
D = j_data[laser_scan]['featurediameter'] # Durchmesser der Linien. Zahl muss ungerade sein! Großzügig abschätzen.
          # diameter = (D,L) ist die Ausdehnung der Objekte in y- und x-Richtung. Tupel-Reihenfolge: (y,x).
diameter = (L,D);
m = j_data[laser_scan]['minmass'] #150 minmass (der kleinste mass-Wert (Maß für Helligkeit), den ein Objekt haben kann).
sep = (D+1,D+1) # separation: der Mindestabstand zwischen den Objekten in y- und x-Richtung. Tupel-Reihenfolge: (y,x).
N = j_data[laser_scan]['serieslen']
first_frame = 0 # Python fängt bei 0 an zu zählen, d. h. frames[0] ist Bild-Nr. 1.
last_frame = first_frame + N-1  # Um das Geschwindigkeitsfeld der Mikropartikel zu berechnen, betrachten wir frames[first_frame] bis 
                # frames[last_frame].
#f = tp.locate(frames[0], diameter=(D,L), minmass=m, separation=sep)
f = tp.locate(frames[first_frame], diameter=(D,L), separation=sep, minmass=m);

# Erstellung eines Balkendiagramms: Anzahl vs mass.  
fig, ax = plt.subplots(figsize=(10,6))
f1 = f.dropna(subset=['mass']) # Herausfiltern der NaN-Werte.
ax.hist(f1['mass'], bins = 20)
ax.set(xlabel='mass', ylabel='Anzahl')
plt.show() 

# Anzeige von frames[first_frame], wobei alle lokalisierten Objekte mit roten Kreisen markiert sind.
plt.figure(figsize=(14,14))
tp.annotate(f, frames[first_frame]);

# An dieser Stelle kann man überprüfen, wie gut die Parameter L, D, m und sep gewählt wurden.


# Bestimmung der Positionen der Mikropartikel auf den Bildern frames[first_frame] bis frames[last_frame]:
f = tp.batch(frames[first_frame:N], diameter=(D,L), minmass=m, separation=sep)
f = f.dropna(axis=0, how='any')

# Verbinden der Positionen der Mikropartikel von den Bildern frames[first_frame] bis frames[last_frame]:
t = tp.link_df(f, search_range=L, adaptive_stop=D, adaptive_step=0.95, memory=0)




# Herausfiltern von kurzlebigen Trajektorien, die nur über wenige Bilder hinweg auftreten. Sie sind meistens keine echten
# Trajektorien.
t1 = tp.filter_stubs(t, N) # Das zweite Argument gibt die minimale Anzahl an Bildern an, die eine Trajektorie andauern muss. 
                           # Wir betrachten nur Trajektorien, die über die gesamte Anzahl an Bildern N andauern. 

# Anzahl der Partikel in den ungefilterten und gefilterten Daten vergleichen:
print('Before:', t['particle'].nunique())
print('After:', t1['particle'].nunique())

# Anzeige von frames[first_frame] mit den lokalisierten Objekten nach der Filterung. 
plt.figure(figsize=(14,14))
tp.annotate(t1[t1['frame'] == 0], frames[0]);


# Bestimmung der horizontalen Geschwindigkeit (Drift enthalten):
data = pd.DataFrame()
for item in set(t1.particle): # set(t1.particle) ist eine Liste der Partikelnummern. 

    sub = t1[t1.particle==item] # Betrachtung eines bestimmten Partikels. t1[t1.particle==item] ist ein DataFrame, in dem die
                                # zur Partikelnummer item zugehörigen Werte für jede frame-Nr. stehen. 
    y_m_px = (sub.y[first_frame]+sub.y[last_frame])/2 # y-Mittelwert in Pixel. Zur Positionsangabe für die x-Geschw. 
    dx = sub.x[last_frame]-sub.x[first_frame] # Zurückgelegte Strecke in x-Richtung.
    v_x = dx*(frame_rate/(N-1))*spatial_resolution_mm_to_pix # x-Geschw. in mm/s (Bildrate = 60 fps).
   
    # DataFrame mit allen benötigten Werten erstellen (y-Mittelwert und x-Geschw. auf 2 Dezimalstellen gerundet):
    data = data.append([{  
                        'y_m_px': np.around(y_m_px, decimals=2),                           
                        'dx': dx, 
                        'v_x': np.around(v_x, decimals=2),
                        'particle': item # Partikel-Nr.
                        }])

data # Zur Überprüfung.



# Erstellung des Geschwindigkeitsprofils (Drift enthalten):
fig, ax = plt.subplots()
font = {'size': 25}
v_x_max = data['v_x'].max() # Benötigen das Maximum der horizontalen Geschwindigkeit zur Achsenskalierung beim Geschwindigkeits-
                            # profil. 
y_m_px_max = data[data.v_x == v_x_max].y_m_px.values[0] # Position des Vektors mit der maximalen horizontalen Geschwindigkeit.
ax.scatter(data.y_m_px*spatial_resolution_mm_to_pix, data.v_x, c='black')
ax.set_xlim((0,frames[0].shape[0]*spatial_resolution_mm_to_pix))
ax.set_ylim((0,v_x_max+0.2))
ax.tick_params(labelsize=25)
ax.set_xlabel('$y$ [mm]', fontdict=font)
ax.set_ylabel('Horizontale Geschwindigkeit [mm/s]', fontdict=font)
plt.show()
print('Maximum der horizontalen Geschwindigkeit:', v_x_max, 'mm/s', 'bei y =',  np.round(y_m_px_max*spatial_resolution_mm_to_pix,2), 'mm')



drift = tp.compute_drift(t1, N) # Berechnung des Partikeldrifts. Das zweite Argument ist ein Glättungsparameter, der die
                                # Anzahl der Bilder angibt, die zur Glättung des Drifts benutzt werden sollen. Es sollen alle N
                                # Bilder benutzt werden. 
t1_corrected = tp.subtract_drift(t1.copy(), drift) # Abzug des Partikeldrifts.
print('Von trackpy berechneter Drift in mm/s:')
drift.x.values[first_frame]*(frame_rate/(N-1))*spatial_resolution_mm_to_pix



# Bestimmung der horizontalen Geschwindigkeit mit Abzug des Drifts:
data_corrected = pd.DataFrame()
for item in set(t1_corrected.particle): # set(t1_corrected.particle) ist eine Liste der Partikelnummern. 
    sub = t1_corrected[t1_corrected.particle==item] # Betrachtung eines bestimmten Partikels. 
                                                    # t1_corrected[t1_corrected.particle==item] ist ein DataFrame, in dem die
                                                    # zur Partikelnummer item zugehörigen Werte für jede frame-Nr. stehen. 
    y_m_px = (sub.y[first_frame]+sub.y[last_frame])/2 # y-Mittelwert in Pixel. Zur Positionsangabe für die x-Geschw. 
    dx = sub.x[last_frame]-sub.x[first_frame] # Zurückgelegte Strecke in x-Richtung.
    v_x = dx*(frame_rate/(N-1))*spatial_resolution_mm_to_pix # x-Geschw. in mm/s (Bildrate = 60 fps).
   
    # DataFrame mit allen benötigten Werten erstellen (y-Mittelwert und x-Geschw. auf 2 Dezimalstellen gerundet):
    data_corrected = data_corrected.append([{
                                            'y_m_px': np.around(y_m_px, decimals=2),                         
                                            'dx': dx, 
                                            'v_x': np.around(v_x, decimals=2),
                                            'particle': item # Partikel-Nr.
                                            }])

data_corrected # Zur Überprüfung.



# Erstellung des Geschwindigkeitsprofils mit Abzug des Drifts:
fig, ax = plt.subplots()
font = {'size': 25}
v_x_max = data_corrected['v_x'].max()  
y_m_px_max = data_corrected[data_corrected.v_x == v_x_max].y_m_px.values[0]
ax.scatter(data_corrected.y_m_px*spatial_resolution_mm_to_pix, data_corrected.v_x, c='black')
ax.set_xlim((0,frames[0].shape[0]*spatial_resolution_mm_to_pix))
ax.set_ylim((0,v_x_max+0.2))
ax.tick_params(labelsize=25)
ax.set_xlabel('$y$ [mm]', fontdict=font)
ax.set_ylabel('Horizontale Geschwindigkeit [mm/s]', fontdict=font)
plt.show()
print('Maximum der horizontalen Geschwindigkeit:', v_x_max, 'mm/s', 'bei y =',  np.round(y_m_px_max*spatial_resolution_mm_to_pix,2), 'mm')



# Manuelle Drift-Berechnung:
dc = data[data['v_x'] > 0] # Partikel, die sich in negative x-Richtung bewegen, werden rausgenommen. Dabei handelt es sich
                           # meistens um keine echten Trajektorien, sondern um von trackpy falsch verbundene Partikel. 
dc_l = dc[dc['y_m_px'] < 300] # Driftende Partikel oberhalb des Laserstrahls (linker Bereich auf der y-Achse).
dc_r = dc[dc['y_m_px'] > 700] # Driftende Partikel unterhalb des Laserstrahls (rechter Bereich auf der y-Achse).
drift_left = dc_l['v_x'].mean()
drift_right = dc_r['v_x'].mean()
# Berechnung der mittleren Driftgeschwindigkeit in x-Richtung: Bei Videos mit sehr wenigen Partikeln kann es sein, dass sich
# oberhalb oder unterhalb des Laserstrahls gar keine driftenden Partikel befinden.
if np.isnan(drift_left) == False and np.isnan(drift_right) == True:
    meandrift = drift_left
elif np.isnan(drift_left) == True and np.isnan(drift_right) == False:
    meandrift = drift_right
elif np.isnan(drift_left) == False and np.isnan(drift_right) == False:
    meandrift = 0
elif np.isnan(drift_left) == True and np.isnan(drift_right) == True:
    meandrift = (drift_left + drift_right)/2

print('Manuell berechneter Drift in mm/s:')
meandrift



# DataFrame mit allen benötigten Werten erstellen:
data_corrected_manually=pd.DataFrame({ 
                                'y_m_px': data.y_m_px, 
                                'v_x': np.around((data.v_x - meandrift), decimals=2), # x-Geschw. abzüglich des Drifts.
                                'particle': data.particle # Partikel-Nr.
                                })
data_corrected_manually # Zur Überprüfung.



# Erstellung des Geschwindigkeitsprofils mit Abzug des manuell berechneten Drifts:
fig, ax = plt.subplots()
font = {'size': 25}
v_x_max = data_corrected_manually['v_x'].max()  
y_m_px_max = data_corrected_manually[data_corrected_manually.v_x == v_x_max].y_m_px.values[0]
ax.scatter(data_corrected_manually.y_m_px*spatial_resolution_mm_to_pix, data_corrected_manually.v_x, c='black')
ax.set_xlim((0,frames[0].shape[0]*spatial_resolution_mm_to_pix))
ax.set_ylim((0,v_x_max+0.2))
ax.tick_params(labelsize=25)
ax.set_xlabel('$y$ [mm]', fontdict=font)
ax.set_ylabel('Horizontale Geschwindigkeit [mm/s]', fontdict=font)
plt.show()
print('Maximum der horizontalen Geschwindigkeit:', v_x_max, 'mm/s', 'bei y =',  np.round(y_m_px_max*spatial_resolution_mm_to_pix,2), 'mm')

endzeit = time.time()
print('Gesamtzeit: {:5.3f}s'.format(endzeit-startzeit))