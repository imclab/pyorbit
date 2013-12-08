pyorbit
=======

Simple orbital dynamics simulator

```python
  # Set up system, will be Earth-centered inertial frame (ECI)
  sys = Orbits()
  # timestep: every hour
  sys.dt = 60*60
  # amount of time to simulate: 28 days
  sys.tfinal = 60 * 60 * 24 * 28
  
  # Add the Earth and Moon (m-k-s units for initial position and velocity)
  sys.add_body(5.97219e24, [0, 0], [0, 0])  # Earth
  sys.add_body(7.3477e22, [405503000, 0], [0, 964.])  # Moon

  # Place a craft in a geostationary orbit (~35,000 km)
  sys.add_craft([-35786000, 0], [0, -3340])
  # Place small craft into a lower orbit
  sys.add_craft([-4140000, 0], [0, -7700])

  # apply some impulse (kind of unstable at the moment...)
  # boost first sattellite towards moon
  sys.impulses[0][350] = 16000
  # boost 2nd to geostationary
  sys.impulses[1][75] = 3500

  # solve the system
  sys.solve_crafts()
  
  # create figure of results
  sys.plot()
  pylab.show()
```

![Alt text](http://i.imgur.com/4M0nTw4.png "Screenshot")
