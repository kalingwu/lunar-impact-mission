function full_project_2
close all
% Satellite's Orbit Calculation
GM=3.986012*10^5; 
theta=7*pi/180;
rp=250+6378.145;
ra=35950+6378.145;
omega=178*pi/180;
ohm=pi;
time=24*60*60;

X=rp*(cos(omega)*cos(ohm)-sin(omega)*sin(ohm)*cos(theta));
Y=rp*(cos(omega)*sin(ohm)+sin(omega)*cos(ohm)*cos(theta));
Z=rp*sin(omega)*sin(theta);

Vx=((2*GM*ra/(rp*(ra+rp)))^0.5)*-(sin(omega)*cos(ohm)+cos(omega)*sin(ohm)*cos(theta));
Vy=((2*GM*ra/(rp*(ra+rp)))^0.5)*(cos(omega)*cos(ohm)*cos(theta)-sin(omega)*sin(ohm));
Vz=((2*GM*ra/(rp*(ra+rp)))^0.5)*(cos(omega)*sin(theta));

w0 = [X,Y,Z,Vx,Vy,Vz]; % Initial conditions
w0pergeeoriginal = w0;
options = odeset('RelTol',0.000000001);

[t_values,w_values] = ode45(@sateom,[0,time],w0,options); 
figure
plot3(w_values(:,1),w_values(:,2),w_values(:,3));

%Moon's Orbit Calculation
%2011 Mar 19 19:00:00.0     11 48 56.844     -  4 15  2.71      356575.001
%2011 Mar 25 07:00:00.0     17 13 41.530     - 23 41 59.10      379613.852
%2011 Apr 02 09:00:00.0     23 47 21.481     +  4  4 52.97      406656.171
GM=3.986012*10^5;
pp=356575.001;
pa=406656.171;
pb=379613.852;
delta=-(4+15/60+2.71/3600)*2*pi/360;
alpha=(11+48/60+56.844/3600)*2*pi/24;
delta_1=-(23+41/60+59.10/3600)*2*pi/360;
alpha_1=(17+13/60+41.530/3600)*2*pi/24;
time=30*24*60*60+8*60*60;

Xmoon=pp*cos(delta)*cos(alpha);
Ymoon=pp*cos(delta)*sin(alpha);
Zmoon=pp*sin(delta);
rp=[Xmoon,Ymoon,Zmoon];

A=pb*cos(delta_1)*cos(alpha_1);
B=pb*cos(delta_1)*sin(alpha_1);
C=pb*sin(delta_1);
rb=[A,B,C];

m=cross(rp,rb);
m = m/sqrt(dot(m,m));
V=(2*GM*pa/(pp*(pa+pp)))^0.5*cross(m,rp)/pp;

wmoon0 = [[Xmoon,Ymoon,Zmoon],V];%initial condition
iccheck = dot([Xmoon,Ymoon,Zmoon],V);

%Special points where the moon's orbit crosses the plane of the satellite
%orbit.
options = odeset('RelTol',0.000000001, 'Event', @detect_impact_point);

[t_moon_values,w_moon_values,tevent,wevent,indexevent] = ode45(@sateom,[0,time],wmoon0,options); 
tevent_acc0=tevent
wevent
indexevent;

hold on
plot3(w_moon_values(:,1),w_moon_values(:,2),w_moon_values(:,3));
rmoon_1=wevent(1,1:3);
hold on
plot3(rmoon_1(1),rmoon_1(2),rmoon_1(3),'ro','MarkerSize',11,'MarkerFaceColor','r');
rmoon_2=wevent(2,1:3);
hold on
plot3(rmoon_2(1),rmoon_2(2),rmoon_2(3),'ro','MarkerSize',11,'MarkerFaceColor','b');


%Trajectory of the launching satellite to the perigee
%[test_distance,test_time]=mindisttomoon(0.67,w0(1:3),w0(4:6),rmoon_2)

    for i =1:41
    delta_V(i)=0.4+0.6*(i-1)/41;
    [dist(i),traveltime(i)]=mindisttomoon(delta_V(i),w0pergeeoriginal(1:3),w0pergeeoriginal(4:6),rmoon_2);
    end
    
    figure
    plot(dist,delta_V);
    [distmin,index]=min(dist);
    critical_delta_V_perigee=delta_V(index);
    figure
    plot(delta_V,traveltime);
    time_to_hit_from_perigee=traveltime(index);

%The critical value of delta_V is 0.6927km/s. The time taken to reach the
%impact point is  2.3349e+05sec.
    
%Trajectory of the launching satellite to the apogee
GM=3.986012*10^5; 
theta=7*pi/180;
rp=250+6378.145;
ra=35950+6378.145;
omega=178*pi/180;
ohm=pi;
time=24*60*60;

X=rp*(cos(omega)*cos(ohm)-sin(omega)*sin(ohm)*cos(theta));
Y=rp*(cos(omega)*sin(ohm)+sin(omega)*cos(ohm)*cos(theta));
Z=rp*sin(omega)*sin(theta);

Vx=((2*GM*ra/(rp*(ra+rp)))^0.5)*-(sin(omega)*cos(ohm)+cos(omega)*sin(ohm)*cos(theta));
Vy=((2*GM*ra/(rp*(ra+rp)))^0.5)*(cos(omega)*cos(ohm)*cos(theta)-sin(omega)*sin(ohm));
Vz=((2*GM*ra/(rp*(ra+rp)))^0.5)*(cos(omega)*sin(theta));

w0 = [X,Y,Z,Vx,Vy,Vz]; % Initial conditions of apogee of the satellite

rapogee=-(ra)*w0(1:3)/rp;
Va=-w0(4:6)*(rp)/ra;
figure
plot3(w_values(:,1),w_values(:,2),w_values(:,3));
hold on
plot3(w_moon_values(:,1),w_moon_values(:,2),w_moon_values(:,3));
hold on
plot3(rmoon_1(1),rmoon_1(2),rmoon_1(3),'ro','MarkerSize',11,'MarkerFaceColor','r');
hold on
plot3(rmoon_2(1),rmoon_2(2),rmoon_2(3),'ro','MarkerSize',11,'MarkerFaceColor','b');
    
    for i =1:41
    delta_V(i)=1+2*(i-1)/41;
    [dist(i),traveltime(i)]=mindisttomoon(delta_V(i),rapogee,Va,rmoon_1);
    end
    
    figure
    plot(dist,delta_V);
    [distmin,index]=min(dist);
    critical_delta_V_apogee=delta_V(index);
    figure
    plot(delta_V,traveltime);
    time_to_hit_from_apogee=traveltime(index);
    
%More accurate calculation of the impact velocity with the moon.
time=30*24*3600;

options = odeset('RelTol',0.000000001);
% [t_test_values,w_test_values] = ode45(@sat_moon_eom,[0,time],[wmoon0,w0pergeeoriginal],options);



options = odeset('RelTol',0.000000001);
t1=tevent_acc0(2)
t2=time_to_hit_from_perigee
[t_moon_accurate,w_moon_accurate] = ode45(@sateom,[0,t1-t2],wmoon0,options);

w_accurate0(1:3)=w_moon_accurate(end,1:3);
w_accurate0(4:6)=w_moon_accurate(end,4:6);
w_accurate0(7:9)=w0pergeeoriginal(1:3);
w_accurate0(10:12)=w0pergeeoriginal(4:6)+0.6927*w0pergeeoriginal(4:6)/sqrt(dot(w0pergeeoriginal(4:6),w0pergeeoriginal(4:6)));

time=30*24*3600;
options = odeset('RelTol',0.001,'Event', @impact);
[t_impact,ws,tevent_impact,wevent_impact,indexevent] = ode45(@sat_moon_eom,[0,time],w_accurate0,options);
tevent_impact;
wevent_impact;

impact_velocity=norm(wevent_impact(end,4:6)-wevent_impact(end,10:12))

for i=1:length(t_impact)
    clf %Clear the frame
    plot3(ws(:,1),ws(:,2),ws(:,3)) %Plot the entire trajectory
    hold on
    plot3(ws(i,1),ws(i,2),ws(i,3),'ro','MarkerSize',10,'MarkerFaceColor','r')
    plot3(ws(:,7),ws(:,8),ws(:,9)) %Plot the entire trajectory
    hold on
    plot3(ws(i,7),ws(i,8),ws(i,9),'ro','MarkerSize',10,'MarkerFaceColor','b')
    pause(0.05)
end

    
    %Function for the equation of motion of the satellite and the moon.
    function dwdt = sateom(t,w)
    x=w(1); y=w(2); z=w(3);
    vx=w(4); vy=w(5); vz=w(6);
    r = sqrt(x^2+y^2+z^2);
    dxdt=vx; dydt=vy; dzdt=vz;
    dvxdt=-GM*x/r^3;
    dvydt=-GM*y/r^3;
    dvzdt=-GM*z/r^3;
    dwdt = [dxdt;dydt;dzdt;dvxdt;dvydt;dvzdt]; 
    end

    for i =1:length(t_values);
    r = sqrt(w_values(i,1)^2 + w_values(i,2)^2+w_values(i,3)^2); 
    vmag = sqrt(w_values(i,4)^2 + w_values(i,5)^2+w_values(i,6)^2);
    
    energy(i) = -GM/r + vmag^2/2; 
    
    a=[w_values(i,1);w_values(i,2);w_values(i,3)];
    b=[w_values(i,4);w_values(i,5);w_values(i,6)];
    angularm(i,1:3) = cross(a,b);
   
    end

    figure
    plot(t_values,energy);
    
    figure
    plot(t_values,angularm);
    
    for i =1:length(t_moon_values);
    R = sqrt(w_moon_values(i,1)^2 + w_moon_values(i,2)^2+w_moon_values(i,3)^2);
    vmag = sqrt(w_moon_values(i,4)^2 + w_moon_values(i,5)^2+w_moon_values(i,6)^2);
    distance(i)=R;
    end
    
    figure
    plot(t_moon_values,distance);
    %From the graph, the time required for one complete orbit of the moon
    %the moon is 2.346*10^6 seconds, which is 27.153 days.
    %Further more, the distance of the moon from the earth at its perigee
    %is 3.566*10^5km and that at its apogee is 4.067*10^5km.
    
    %Function to detect the impact point when the moon is in the plane of
    %the satellite's orbit.
    function [event_val,stopthecalc,direction]=detect_impact_point(t_moon_values,w_moon_values)
        %formula for unit vector normal to satellite orbit plane
        n=[sin(omega)*sin(theta),-cos(omega)*sin(theta),cos(theta)];
        %Position vector of moon (this is auumes w_moon_values(1)=x, 
        %w_moon_values(2)=y, w_moon_values(3)=z)
        rr=w_moon_values(1:3);
        %Detect when r.n=0
        event_val=dot(rr,n);
        stopthecalc=0;
        direction=0;
    end

    %function to calculate min distance to moon for satellite trajectory.
    function [d,time_to_reach_min_point]=mindisttomoon(delta_V,r0,v0,r_impact)
    time=8*24*60*60;
    w0(1:3)=r0;
    w0(4:6)=v0+delta_V*v0/sqrt(dot(v0,v0));
    rtol=0.00000001
    options=odeset('Reltol',rtol,'Events',@min_dist);
    [t_vals,w_vals,tevent,wevent,ievent]=ode45(@sateom,[0,time],w0,options);
    wevent
    hold on
    plot3(w_vals(:,1),w_vals(:,2),w_vals(:,3))
    rmin=wevent(1,1:3)-r_impact(1:3)
    d=sqrt(dot(rmin,rmin));
    time_to_reach_min_point=tevent(1);
   
    function [event_val,stopthecalc,direction]=min_dist(t_vals,w_values)
        r=transpose(w_values(1:3));
        v=transpose(w_values(4:6));
        event_val=dot((r-r_impact),v);
        stopthecalc=1;
        direction=0;
    end
    
    end

    %Function for the equation of motion of the satellite account for the
    %moon's gravity.
    function dwdt = sat_moon_eom(t,w)
    xm=w(1); ym=w(2); zm=w(3);
    vxm=w(4); vym=w(5); vzm=w(6);
    
    xs=w(7); ys=w(8); zs=w(9);
    vxs=w(10); vys=w(11); vzs=w(12);
    
    rm = sqrt(xm^2+ym^2+zm^2);
    rs = sqrt(xs^2+ys^2+zs^2);
    rms=sqrt((xs-xm)^2+(ys-ym)^2+(zs-zm)^2);
    dxmdt=vxm; dymdt=vym; dzmdt=vzm;
    dxsdt=vxs; dysdt=vys; dzsdt=vzs;
    dvxmdt=-GM*xm/rm^3;
    dvymdt=-GM*ym/rm^3;
    dvzmdt=-GM*zm/rm^3;
    dvxsdt=-GM*xs/rs^3-0.012298*GM*(xs-xm)/rms^3;
    dvysdt=-GM*ys/rs^3-0.012298*GM*(ys-ym)/rms^3;
    dvzsdt=-GM*zs/rs^3-0.012298*GM*(zs-zm)/rms^3;
    dwdt = [dxmdt;dymdt;dzmdt;dvxmdt;dvymdt;dvzmdt;dxsdt;dysdt;dzsdt;dvxsdt;dvysdt;dvzsdt];
    end

    function [event_val,stopthecalc,direction]=impact(t,w)
    rs=w(7:9)
    rm=w(1:3)
    ds=sqrt(dot(rs-rm,rs-rm));
    event_val=ds-1737.4;    
    stopthecalc=1;
    direction=0;
    end

end