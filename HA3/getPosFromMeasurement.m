function [x, y] = getPosFromMeasurement(y1, y2, s1, s2)
  %GETPOSFROMMEASUREMENT computes the intersection point
  %(transformed 2D measurement in Cartesian coordinate
  %system) given two sensor locations and two bearing
  %measurements, one from each sensor.
  %INPUT: y1: bearing measurement from sensor 1
  %       y2: bearing measurement from sensor 2
  %       s1: location of sensor 1 in 2D Cartesian
  %       s2: location of sensor 2 in 2D Cartesian
  %OUTPUT: x: coordinate of intersection point on x axis
  %        y: coordinate of intersection point on y axis
  %This problem can be formulated as solving a set of two
  %linear equations with two unknowns. Specifically, one
  %would like to obtain (x,y) by solving
  %(y - s1(2))=(x - s1(1))tan(y1) and (y - s2(2))=(x - s2(1))tan(y2).
  x = (s2(2) - s1(2) + tan(y1)*s1(1) - tan(y2)*s2(1))/(tan(y1) - tan(y2));
  y = s1(2) + tan(y1)*(x - s1(1));
end