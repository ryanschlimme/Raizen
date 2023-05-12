# 11 May 2023
# Ryan Schlimme

# Starting to explore ray tracing to model deflection of laser rays through wavefronts of different refractive index
# Goals: develop method to have radially expanding circle in time, then moving circle center in time, concentric circles of unique refractive index, continuous wave of refractive index in time

import numpy as np
import matplotlib.pyplot as plt

def quadratic(x, k):
    return 1/2 * k * x*x

def normalize(vector):
    return vector / np.linalg.norm(vector)

class Ray:
    def __init__(self, origin, direction, power, polarization):
        self.origin = origin
        self.direction = direction
        self.power = power
        self.polarization = polarization

    def termination(self, length):
        return self.origin + self.direction * length
    
    def plot(self, length=None, ax=None, arrowstyle="-", **kwargs):
        if length is None:
            termination = self.intersection
        else:
            termination = self.termination(length)
        if ax is None:
            fig, ax = plt.subplots(1,1)
        ax.annotate("", xytext=self.origin, xy=termination, annotation_clip=False,
                    arrowprops=dict(arrowstyle=arrowstyle, alpha=np.sqrt(self.power), lw=1)|kwargs)
    
    def sphere_intersect(self, center, radius, origin=None):
        if origin is None:
            origin = self.origin
        b = 2 * np.dot(self.direction, origin - center)
        c = np.linalg.norm(origin - center) ** 2 - radius ** 2
        delta = b ** 2 - 4 * c
        if delta > 0:
            d1 = (-b + np.sqrt(delta)) / 2
            d2 = (-b - np.sqrt(delta)) / 2
            d1, d2 = np.sort([d1, d2])
            return d1, d2
        else:
            return None, None
        
    def first_intersection(self, center, radius):
        distance, _ = self.sphere_intersect(center, radius)
        if distance is None:
            self.intersection = None
        return self.termination(distance)
        
    def next_intersection(self, center, radius):
        origin_normal = normalize(self.origin - center)
        shifted_origin = self.origin + 1e-5 * origin_normal
        _, chord_length = self.sphere_intersect(center, radius, origin=shifted_origin)
        return self.termination(chord_length)
    
    def on_sphere(self, center, radius):
        return (np.linalg.norm(self.origin - center) - radius) < 1e-5
    
    def set_intersection(self, center, radius):
        if self.on_sphere(center, radius):
            self.intersection = self.next_intersection(center, radius)
        else:
            self.intersection = self.first_intersection(center, radius)       
    
    def sphere_normal(self, center):
        return normalize(self.intersection - center)
    
    def cosin_angles(self, center, r):
        normal = self.sphere_normal(center)
        c1 = -np.dot(self.direction, normal)
        if c1 < 0:
            normal = -normal
            c1 = -c1
        c2 = np.sqrt(np.abs(1 - r*r*(1 - c1*c1)))
        return c1, c2, normal
    
    def reflected_direction(self, center, r=1):
        c1, _, normal = self.cosin_angles(center, r)
        return self.direction + 2 * c1 * normal

    def refracted_direction(self, center, r):
        c1, c2, normal, = self.cosin_angles(center, r)
        return r*self.direction + (r*c1 - c2) *  normal
    
    def reflection_coefficient(self, center, r):
        normal = self.sphere_normal(center)
        c1, c2, normal = self.cosin_angles(normal, r)
        if self.polarization == "s":
            return np.abs((r*c1 - c2)/(r*c1 + c2))**2
        if self.polarization == "p":
            return np.abs((r*c2 - c1)/(r*c2 + c1))**2
    
    def transmission_coefficient(self, center, r):
        return 1 - self.reflection_coefficient(center, r)
    
    def reflected_power(self, center, r):
        return self.power * self.reflection_coefficient(center, r)
    
    def transmitted_power(self, center, r):
        return self.power * self.transmission_coefficient(center, r)
    
    def spawn_transmitted_ray(self, center, r):
        return Ray(origin=self.intersection, 
                direction=self.refracted_direction(center, r), 
                power=self.transmitted_power(center, r),
                polarization=self.polarization
               )

    def spawn_reflected_ray(self, center, r):
        return Ray(origin=self.intersection, 
                direction=self.reflected_direction(center, r), 
                power=self.reflected_power(center, r),
                polarization=self.polarization
               )
    
    def spawn_rays(self, center, r):
        reflected_ray = self.spawn_reflected_ray(center, r)
        transmitted_ray = self.spawn_transmitted_ray(center, r)
        return reflected_ray, transmitted_ray
    

def trace_rays(ray0, center, radius, n1, n2, N_reflections):
    ray0.set_intersection(center, radius)
    if ray0.intersection is None:
        return None
    rayR1, rayT1 = ray0.spawn_rays(center, n1/n2)
    rayT1.set_intersection(center, radius)
    internal_rays = [rayT1]
    external_rays = [rayR1]
    for i in range(N_reflections):
        internal_ray, external_ray = internal_rays[i].spawn_rays(center, n2/n1)
        internal_ray.set_intersection(center, radius)
        internal_rays.append(internal_ray)
        external_rays.append(external_ray)
    
    external_ray = internal_rays[-1].spawn_transmitted_ray(center, n2/n1)
    external_rays.append(external_ray)
    return external_rays, internal_rays


fig, ax = plt.subplots(1,1, figsize=(5.5, 3))

radius = 1
focal_length  = 6*radius
n1 = 1.0
n2 = 1.47
c = 3e8
beam_power = 50e-3
beam_waist = radius
polarization = "s"
N_reflections = 0

R_lens = 3*radius
N_rays = 10
x_lens = np.linspace(-R_lens, R_lens,  N_rays)

center=[focal_length, -0.5]

P0 = beam_power / len(x_lens)
#P0 = beam_power / np.sum(np.exp(-2*x_lens**2/beam_waist**2))
force = 0
for x0 in x_lens:
    focus = np.array([focal_length, 0])
    #focus = np.array([focal_length, x0])
    origin = np.array([0, x0])
    direction = normalize(focus - origin)
    power = P0
    #power = P0 * np.exp(-2*(x0**2)/beam_waist**2)
    ax.annotate("", xytext=[-radius, x0], xy=[0, x0], annotation_clip=False,
                arrowprops=dict(arrowstyle="-", color="k"))
    #rayp.plot(ax=ax, length=2.3, alpha=np.sqrt(ray0.power/P0),arrowstyle="-")
    ray0 = Ray(origin, direction, power, polarization=polarization) 
    try:
        external_rays, internal_rays = trace_rays(ray0, center, radius, n1, n2, N_reflections)
    except TypeError:
        if ax:
            ray0.plot(ax=ax, length=2.3, alpha=np.sqrt(ray0.power/P0),arrowstyle="-")
        continue       
    force += n2/c * (ray0.power*ray0.direction - sum(ray.power*ray.direction for ray in external_rays))   
    if ax:
        ray0.plot(ax=ax, alpha=np.sqrt(ray0.power/P0), arrowstyle="-")    
        for ray in internal_rays:
            ray.plot(ax=ax,alpha=np.sqrt(ray.power/P0))
        for rayi, ray in enumerate(external_rays):
            if rayi ==1:
                length = 6*radius
                ax.annotate("", xytext=ray.termination(length), xy=np.array([14, ray.termination(length)[1]]), annotation_clip=False,
                arrowprops=dict(arrowstyle="-", color="k"))
                ray.plot(ax=ax, alpha=np.sqrt(ray.power/P0), length=length, arrowstyle="-")
    
if ax:
    zs_circle = np.linspace(center[0]-radius, center[0]+radius, 1000)
    xs_circle1 = center[1] + np.sqrt(radius**2 - (zs_circle-center[0])**2)
    xs_circle2 = center[1] - np.sqrt(radius**2 - (zs_circle-center[0])**2)
    ax.plot(zs_circle, xs_circle1, c="k")
    ax.plot(zs_circle, xs_circle2, c="k")
    #ax.annotate("", xytext=center, xy=center+force*5e10,annotation_clip=False, 
    #            arrowprops=dict(arrowstyle="->", color="crimson"))
ax.set_xlim(focal_length-8*radius, focal_length+7*radius)
ax.set_ylim(-4*radius, 4*radius)
ax.set_aspect(1)
ax.set_axis_off()
plt.tight_layout()
plt.show()