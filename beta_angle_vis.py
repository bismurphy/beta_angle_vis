from skyfield.api import load, EarthSatellite,wgs84
from skyfield.elementslib import osculating_elements_of
import load_tle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider,Button

tle = load_tle.get_tle(53768)
year = 2022
start_day = 180

sat = EarthSatellite(*tle)
eph = load('de421.bsp')
earth,sun = eph['earth'], eph['sun']
ts = load.timescale()

#Plot the data. Points in front of earth show in black, points behind in red (like in finances)
fig,ax = plt.subplots()
earth_circle = plt.Circle((0,0),radius=6371,color='blue')
ax.add_patch(earth_circle)
plt.subplots_adjust(bottom=0.25)

date_slider = Slider(
    ax = plt.axes([0.2, 0.1, 0.6, 0.03]),
    label="Date",
    valmin=0,
    valmax=365,
    valinit=start_day,
    valstep=1
)
prev_year = Button(
    plt.axes([0.2, 0.05, 0.15, 0.04],facecolor='k'),
    "Last year")
next_year = Button(
    plt.axes([0.7, 0.05, 0.15, 0.04],facecolor='k'),
    "Next year")
prev_year.on_clicked(lambda x:year_shift(-1))
next_year.on_clicked(lambda x:year_shift(1))
def year_shift(amount):
    global year
    year += amount
    date_slider.set_val(1)
def date_update(new_date):
    date_of_interest = [year,1,new_date]
    t = ts.utc(*date_of_interest)
    date_slider.valtext.set_text(t.utc_strftime('%Y-%m-%d'))
    earthat = earth.at(t)
    earth_r = earthat.position.km
    earth_v = earthat.velocity.km_per_s
    #To render the sat from the POV of the sun, we need a coordinate system.
    #X axis points away from the sun (or, from the sun to the earth)
    #Y axis points to the ecliptic north pole
    #Z axis completes the orthogonal set
    #And this will be centered on the earth.
    #Take our established positions and dot them with those vectors.

    #Unit vector that points to the earth
    x_hat = earth_r / np.linalg.norm(earth_r)
    #To get the ecliptic north pole, get the angular momentum vector of earth orbit
    earth_h = np.cross(earth_r,earth_v)
    y_hat = earth_h / np.linalg.norm(earth_h)
    #Finally for z. Just x cross y.
    z_hat = np.cross(x_hat,y_hat)

    #We can make a transformation matrix. satpos dotted with each of these is our output.
    transformation_matrix = np.vstack((x_hat,y_hat,z_hat))
    #Add earth's axis
    north_pole = wgs84.latlon(90,0).itrs_xyz.km
    south_pole = wgs84.latlon(-90,0).itrs_xyz.km
    north_pole =np.dot(transformation_matrix,north_pole)
    south_pole =np.dot(transformation_matrix,south_pole)
    earth_axis.set_xdata([north_pole[2],south_pole[2]])
    earth_axis.set_ydata([north_pole[1],south_pole[1]])

    #Now convert sat position into these coordinates. Run for one period.
    sat_state = sat.at(t)
    sat_period_minutes = osculating_elements_of(sat_state).period_in_days*1440
    one_period = ts.utc(*date_of_interest,0,range(int(sat_period_minutes)+5))
    #This gets 92 x values, 92 y values, and 92 z values. I want 92 positions.
    satpos = sat.at(one_period).position.km
    #This zip operation fixes that
    position_list = list(zip(*satpos))
    #Now matrix multiply.
    satpos_in_sunframe = [list(np.dot(transformation_matrix,pos)) for pos in position_list]
    new_fronts = []
    new_backs = []
    for pos in satpos_in_sunframe:
        if pos[0] > 0:
            new_backs.append(pos[2:0:-1])
        else:
            new_fronts.append(pos[2:0:-1])
    back_scatter.set_offsets(new_backs)
    front_scatter.set_offsets(new_fronts)
    fig.canvas.draw_idle()
earth_axis, = ax.plot([],[],linewidth=2)
back_scatter = ax.scatter([],[],color='red')
front_scatter = ax.scatter([],[],color='black')
date_slider.on_changed(date_update)
ax.axis('equal')
date_update(start_day)
plt.show()
