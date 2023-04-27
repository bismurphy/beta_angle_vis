from skyfield.api import load, EarthSatellite,wgs84
from skyfield.elementslib import osculating_elements_of
import load_tle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider,Button

#53768 is beavercube. NPP is 37849. GOES-18 is 51850. TESS is 43435
tle = load_tle.get_tle(53768)
ts = load.timescale()
#Automatically load up today.
this_year,today_day_of_year=[int(x) for x in (ts.now().utc_strftime("%Y,%-j")).split(",")]
#If you'd like to run with a different day, that option is here
year = this_year
start_day = today_day_of_year

sat = EarthSatellite(*tle)
eph = load('de421.bsp')
earth,sun = eph['earth'], eph['sun']

#Plot the data. Points in front of earth show in black, points behind in red (like in finances)
fig,ax = plt.subplots()
required_plot_size = osculating_elements_of(sat.at(ts.utc(year,1,start_day))).apoapsis_distance.km
required_plot_size *= 1.2
ax.set_xlim([-required_plot_size,required_plot_size])
ax.set_ylim([-required_plot_size,required_plot_size])

ax.set_aspect('equal')
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
    earth_axis_upper.set_xdata([north_pole[2],0])
    earth_axis_lower.set_xdata([south_pole[2],0])
    earth_axis_upper.set_ydata([north_pole[1],0])
    earth_axis_lower.set_ydata([south_pole[1],0])
    #if north pole is forward
    if north_pole[0] < 0:
        earth_axis_upper.set_zorder(0.7)
        earth_axis_lower.set_zorder(0.3)
    else:
        earth_axis_upper.set_zorder(0.3)
        earth_axis_lower.set_zorder(0.7)
    #Generate an equator. At spring/fall, it's a diagonal line. At summer/winter, it's an ellipse that goes up and down by a bit
    #This is how much the axis tilts to the left from vertical (convention is counter clockwise angle is positive)
    visible_axial_tilt = -np.arctan2(north_pole[2],north_pole[1])
    #This is how much the axis tilts away from the camera (on the north pole end)
    axial_tilt_away = np.arctan2(*north_pole[:2])
    #Use these to make a transformation of a circle. start with a set of points for the circle.
    front_equator_points = []
    back_equator_points = []
    equator_height = np.sin(axial_tilt_away)
    for angle in np.linspace(0,np.pi*2,90):
        #Start by drawing a circle around the earth.
        x_value = 6371*np.cos(angle)
        y_value = 6371*np.sin(angle)
        #Squish it down. We calculated height of equator already, so multiply by that.
        y_value *= equator_height
        #Rotate by the visible axis tilt
        x_rot = x_value*np.cos(visible_axial_tilt)-y_value*np.sin(visible_axial_tilt)
        y_rot = x_value*np.sin(visible_axial_tilt)+y_value*np.cos(visible_axial_tilt)
        #We have a fun little thing with the fact that equator height can be negative, turns out the front is always drawn by the first half
        if (angle < np.pi):
            front_equator_points.append([x_rot,y_rot])
        else:
            back_equator_points.append([x_rot,y_rot])
    front_equator.set_xdata(list(zip(*front_equator_points))[0])
    front_equator.set_ydata(list(zip(*front_equator_points))[1])
    back_equator.set_xdata(list(zip(*back_equator_points))[0])
    back_equator.set_ydata(list(zip(*back_equator_points))[1])
    #Earth drawing complete.
    #Now convert sat position into these coordinates. Run for one period.
    sat_state = sat.at(t)
    sat_period_minutes = osculating_elements_of(sat_state).period_in_days*1440
    one_period = ts.utc(*date_of_interest,0,range(int(sat_period_minutes)+1))
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
#Reminder: Z order of 0 is in the back, 1 is in the front! Required order, back to front:
#back scatter (for sat trail), back equator, back axis, earth, front axis, front equator, front scatter
earth_circle = plt.Circle((0,0),radius=6371,color='lightblue',zorder=0.5,alpha=0.8)
ax.add_patch(earth_circle)
#Plot the axis, with a marker on the north pole. Only mark point 0, which is defined on pole.
earth_axis_upper, = ax.plot([],[],linewidth=2,color='orange',marker='.',markevery=[0],markerfacecolor='yellow')
earth_axis_lower, = ax.plot([],[],linewidth=2,color='orange',marker='.',markevery=[0],markerfacecolor='yellow')
front_equator, = ax.plot([],[],linewidth=2,color='green',zorder=0.8)
back_equator, = ax.plot([],[],linewidth=2,color='green',zorder=0.2)
back_scatter = ax.scatter([],[],color='black',zorder=0.1)
front_scatter = ax.scatter([],[],color='red',zorder=0.9)
date_update(start_day)

date_slider.on_changed(date_update)
plt.show()
