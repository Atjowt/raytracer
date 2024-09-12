#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

static const double TAU = 6.283185307179586;
static const double PI = 0.5 * TAU;

double lerp(double a, double b, double t) {
	return a + t * (b - a);
}

double frand(void) {
	return rand() / (double)RAND_MAX;
}

typedef struct {
	double x, y, z;
} vec3;

vec3 vec3_add(vec3 u, vec3 v) {
	return (vec3) { u.x + v.x, u.y + v.y, u.z + v.z };
}

vec3 vec3_sub(vec3 u, vec3 v) {
	return (vec3) { u.x - v.x, u.y - v.y, u.z - v.z };
}

vec3 vec3_scale(vec3 v, double s) {
	return (vec3) { v.x * s, v.y * s, v.z * s };
}

vec3 vec3_prod(vec3 u, vec3 v) {
	return (vec3) { u.x * v.x, u.y * v.y, u.z * v.z };
}

vec3 vec3_lerp(vec3 a, vec3 b, double t) {
	return vec3_add(a, vec3_scale(vec3_sub(b, a), t));
}

double vec3_dot(vec3 u, vec3 v) {
	return u.x * v.x + u.y * v.y + u.z * v.z;
}

double vec3_lensqr(vec3 v) {
	return vec3_dot(v, v);
}

double vec3_len(vec3 v) {
	return sqrt(vec3_lensqr(v));
}

double vec3_distsqr(vec3 a, vec3 b) {
	return vec3_lensqr(vec3_sub(b, a));
}

double vec3_dist(vec3 a, vec3 b) {
	return sqrt(vec3_distsqr(a, b));
}

vec3 vec3_normalize(vec3 v) {
	return vec3_scale(v, 1.0 / vec3_len(v));
}

vec3 vec3_from_angles(double theta, double phi) {
	double sin_phi = sin(phi);
	double cos_phi = cos(phi);
	double sin_theta = sin(theta);
	double cos_theta = cos(theta);
	return (vec3) {
		cos_phi * cos_theta,
		cos_phi * sin_theta,
		sin_phi,
	};
}

vec3 vec3_random(void) {
	double theta = TAU * frand();
	double phi = 0.5 * TAU * frand();
	return vec3_from_angles(theta, phi);
}

vec3 vec3_reflect(vec3 v, vec3 n) {
	return vec3_sub(v, vec3_scale(n, 2.0 * vec3_dot(v, n)));
}

vec3 vec3_refract(vec3 uv, vec3 n, double etai_over_etat) {
    double cos_theta = fmin(-vec3_dot(uv, n), 1.0);
    vec3 r_out_perp = vec3_scale(vec3_add(uv, vec3_scale(n, cos_theta)), etai_over_etat);
    vec3 r_out_parallel = vec3_scale(n, -sqrt(fabs(1.0 - vec3_lensqr(r_out_perp))));
    return vec3_add(r_out_perp, r_out_parallel);
}

static const vec3 UP = { 0.0, 0.0, 1.0 };
static const vec3 RIGHT = { 0.0, 0.0, -1.0 };
static const vec3 FORWARD = { -1.0, 0.0, 0.0 };

static const vec3 RED = { 1.0, 0.0, 0.0 };
static const vec3 GREEN = { 0.0, 1.0, 0.0 };
static const vec3 BLUE = { 0.0, 0.0, 1.0 };
static const vec3 YELLOW = { 0.0, 0.0, 1.0 };
static const vec3 BLACK = { 0.0, 0.0, 0.0 };
static const vec3 WHITE = { 1.0, 1.0, 1.0 };

typedef struct {
	vec3 origin, direction;
} Ray;

Ray ray_between(vec3 a, vec3 b) {
	return (Ray) {
		.origin = a,
		.direction = vec3_normalize(vec3_sub(b, a)),
	};
}

vec3 ray_at(Ray ray, double t) {
	return vec3_add(ray.origin, vec3_scale(ray.direction, t));
}

typedef struct {
	vec3 center;
	double radius;
} Sphere;

typedef struct {
	vec3 point;
	vec3 normal;
} Plane;

typedef struct {
	enum {
		SHAPE_PLANE,
		SHAPE_SPHERE,
	} type;
	union {
		Plane plane;
		Sphere sphere;
	} value;
} Shape;

typedef struct {
	vec3 color;
	double roughness;
	double transparent;
	double metallic;
	double specular;
	double emission;
	double transmission;
	double ior;
} Material;

typedef struct {
	Shape shape;
	Material material;
} Object;

typedef struct {
	bool outside;
	double t;
	vec3 point;
	vec3 normal;
	const Object* object;
} Hit;

typedef struct {
	Object* objects;
	int count;
} World;

bool ray_hit_plane(Ray ray, Plane plane, double t_min, double t_max, Hit* hit) {
	double n_dot_d = vec3_dot(plane.normal, ray.direction);
	if (n_dot_d == 0.0) {
		return false;
	}
	double n_dot_p = vec3_dot(plane.normal, plane.point);
	double n_dot_o = vec3_dot(plane.normal, ray.origin);
	double t = (n_dot_p - n_dot_o) / n_dot_d;
	if (t < t_min || t > t_max) {
		return false;
	}
	if (hit != NULL) {
		hit->t = t;
		hit->point = ray_at(ray, t);
		hit->outside = (n_dot_d < 0.0);
		hit->normal = hit->outside ? plane.normal : vec3_scale(plane.normal, -1.0);
	}
	return true;
}

bool ray_hit_sphere(Ray ray, Sphere sphere, double t_min, double t_max, Hit* hit) {
	/*x^2 + y^2 + z^2 = r^2*/
	/*x = a+t*dx*/
	/*y = b+t*dy*/
	/*z = c+t*dz*/
	/*a^2 + 2*a*t*dx + t^2*dx^2 + b^2 + 2*b*t*dy + t^2*dy^2 + c^2 + 2*c*t*dz + t^2*dz^2 = r^2*/
	/*(dx^2+dy^2+dz^2)*t^2 + (2*a*dx + 2*b*dy + 2*c*dz)*t + a^2+b^2+c^2-r^2 = 0*/
	/*solve for t*/
	double a = ray.origin.x - sphere.center.x;
	double b = ray.origin.y - sphere.center.y;
	double c = ray.origin.z - sphere.center.z;
	double dx = ray.direction.x;
	double dy = ray.direction.y;
	double dz = ray.direction.z;
	double r = sphere.radius;
	double c2 = dx*dx + dy*dy + dz*dz;
	double c1 = 2*(a*dx + b*dy + c*dz);
	double c0 = a*a + b*b +c*c -r*r;
	double p = c1 / c2;
	double q = c0 / c2;
	double discriminant = p*p/4.0 - q;
	if (discriminant < 0.0) {
		return false;
	}
	double offset = sqrt(discriminant);
	double midpoint = -p/2.0;
	double t1 = midpoint - offset;
	double t2 = midpoint + offset;
	double t = t1;
        if (t <= t_min || t_max <= t) {
		t = t2;
		if (t <= t_min || t_max <= t) {
			return false;
		}
	}
	if (hit != NULL) {
		hit->t = t1;
		hit->point = ray_at(ray, t);
		hit->normal = vec3_scale(vec3_sub(hit->point, sphere.center), 1.0 / r);
		hit->outside = (vec3_dot(vec3_scale(ray.direction, -1.0), hit->normal) >= 0.0);
		hit->normal = hit->outside ? hit->normal : vec3_scale(hit->normal, -1.0);
	}
	return true;
}

bool ray_hit_shape(Ray ray, Shape shape, double t_min, double t_max, Hit* hit) {
	switch (shape.type) {
		case SHAPE_PLANE: return ray_hit_plane(ray, shape.value.plane, t_min, t_max, hit);
		case SHAPE_SPHERE: return ray_hit_sphere(ray, shape.value.sphere, t_min, t_max, hit);
	}
	return false;
}

bool ray_hit_object(Ray ray, const Object* restrict object, double t_min, double t_max, Hit* restrict hit) {
	if (ray_hit_shape(ray, object->shape, t_min, t_max, hit)) {
		if (hit != NULL) {
			hit->object = object;
		}
		return true;
	}
	return false;
}

bool ray_hit_world(Ray ray, World world, double t_min, double t_max, Hit* hit) {

	Hit hit_temp;
	bool hit_anything = false;

	for (int i = 0; i < world.count; i++) {
		if (ray_hit_object(ray, &world.objects[i], t_min, t_max, &hit_temp)) {
			hit_anything = true;
			t_max = hit_temp.t;
			if (hit != NULL) {
				*hit = hit_temp;
			}
		}
	}

        return hit_anything;
}

vec3 get_skybox_color(Ray ray) {
	static const vec3 sky_color = { 0.6, 0.7, 0.9 };
	static const vec3 sun_color = { 0.9, 0.8, 0.8 };
	/*static const vec3 sky_color = { 0.1, 0.1, 0.1 };*/
	/*static const vec3 sun_color = { 1.0, 1.0, 1.0 };*/
	vec3 sun_dir = vec3_normalize((vec3) { 0.0, -0.1, -1.0 });
	double sun_fac = fmax(0.0, -vec3_dot(ray.direction, sun_dir));
	return vec3_lerp(sky_color, sun_color, sun_fac);
}

// Schlick's approximation
double get_reflectance(double cos_theta, double ior) {
        double r = (1.0 - ior) / (1.0 + ior);
        double r2 = r * r;
        return r2 + (1.0 - r2) * pow((1.0 - cos_theta), 5);
}

vec3 get_material_color(Material material, Ray ray, Hit hit) {
	/*return hit.outside ? BLUE : RED;*/
	return material.color;
	return vec3_lerp(material.color, WHITE, material.emission);
}

Ray material_bounce_ray(Material material, Ray incoming, Hit hit) {

	Ray outgoing;

	outgoing.origin = hit.point;

	// Diffuse
	/*outgoing.direction = vec3_normalize(vec3_add(hit.normal, vec3_random()));*/
	double theta = TAU * frand();
	double phi = material.roughness * 0.5 * TAU * frand();
	outgoing.direction = vec3_from_angles(theta, phi);
	if (vec3_dot(outgoing.direction, hit.normal) <= 0.0) {
		outgoing.direction = vec3_scale(outgoing.direction, -1.0);
	}

	// Metallic
	vec3 reflected = vec3_reflect(incoming.direction, hit.normal);
        outgoing.direction = vec3_lerp(outgoing.direction, reflected, material.metallic);

	// Transparent
        outgoing.direction = vec3_lerp(outgoing.direction, incoming.direction, material.transparent);

	// Emission
        outgoing.direction = vec3_lerp(outgoing.direction, hit.normal, material.emission);

	// Transmission
	double ior = hit.outside ? (1.0 / material.ior) : material.ior;
	double cos_theta = fmin(-vec3_dot(incoming.direction, hit.normal), 1.0);
	double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
	vec3 refracted;
	if ((ior * sin_theta > 1.0) || (get_reflectance(cos_theta, ior) > frand())) {
		refracted = vec3_reflect(incoming.direction, hit.normal);
	} else {
		refracted = vec3_refract(incoming.direction, hit.normal, ior);
	}
        outgoing.direction = vec3_lerp(outgoing.direction, refracted, material.transmission);

	return outgoing;
}

vec3 trace_ray(Ray ray, World world, int bounces) {
	static const double t_epsilon = 0.001;
	Hit hit;
	vec3 ray_color = WHITE;
	while (true) {
		if (bounces > 0 && ray_hit_world(ray, world, t_epsilon, INFINITY, &hit)) {
			Material material = hit.object->material;
			ray_color = vec3_prod(ray_color, get_material_color(material, ray, hit));
			ray = material_bounce_ray(material, ray, hit);
			bounces--;
		} else {
			return vec3_prod(ray_color, get_skybox_color(ray));
		}
	}
}

double radians_from_degrees(double degrees) {
	return TAU * (degrees / 360.0);
}

double degress_from_radians(double radians) {
	return 360.0 * (radians / TAU);
}

double focal_length_from_fov(double fov_radians) {
	return 0.5 / tan(0.5 * fov_radians);
}

int main(void) {

	static const int image_width = 1920 / 3;
	static const int image_height = 1080 / 3;
	static const double aspect_ratio = image_width / (double)image_height;
	static const double viewport_height = 1.0;
	static const double viewport_width = viewport_height * aspect_ratio;
	static const int max_bounces = 8;
	static const int samples_per_pixel = 16;
	const double fov = radians_from_degrees(60.0);
	const double focal_length = focal_length_from_fov(fov);
	const vec3 eye_pos = { 0.0, -focal_length, 0.0 };

	Object objects[] = {
		(Object) {
			.shape = {
				.type = SHAPE_PLANE,
				.value.plane = {
					.point = { 0.0, 0.0, -0.5 },
					.normal = UP,
				},
			},
			.material = {
				.color = { 0.5, 0.5, 0.5 },
				.roughness = 0.5,
				.metallic = 0.5,
				.transparent = 0.0,
				.emission = 0.0,
				.transmission = 0.0,
				.ior = 1.0,
			},
		},
		(Object) {
			.shape = {
				.type = SHAPE_SPHERE,
				.value.sphere = {
					.center = { -1.1, 1.0, 0.0 },
					.radius = 0.5,
				},
			},
			.material = {
				.color = { 1.0, 0.1, 0.1 },
				.roughness = 0.1,
				.metallic = 1.0,
				.transparent = 0.0,
				.emission = 0.0,
				.transmission = 0.0,
				.ior = 1.0,
			},
		},
		(Object) {
			.shape = {
				.type = SHAPE_SPHERE,
				.value.sphere = {
					.center = { 0.0, 1.0, 0.0 },
					.radius = 0.5,
				},
			},
			.material = {
				.color = { 1.0, 1.0, 1.0 },
				.roughness = 0.0,
				.metallic = 0.0,
				.transparent = 0.0,
				.emission = 0.0,
				.transmission = 1.0,
				.ior = 1.5,
			},
		},
		(Object) {
			.shape = {
				.type = SHAPE_SPHERE,
				.value.sphere = {
					.center = { 0.0, 1.0, 0.0 },
					.radius = 0.4,
				},
			},
			.material = {
				.color = { 1.0, 1.0, 1.0 },
				.roughness = 0.0,
				.metallic = 0.0,
				.transparent = 0.0,
				.emission = 0.0,
				.transmission = 1.0,
				.ior = 1.0 / 1.5,
			},
		},
		(Object) {
			.shape = {
				.type = SHAPE_SPHERE,
				.value.sphere = {
					.center = { 1.1, 1.0, 0.0 },
					.radius = 0.5,
				},
			},
			.material = {
				.color = { 0.1, 0.1, 1.0 },
				.roughness = 0.1,
				.metallic = 1.0,
				.transparent = 0.0,
				.emission = 0.0,
				.transmission = 0.0,
				.ior = 1.0,
			},
		},
	};

	World world = {
		.objects = objects,
		.count = sizeof(objects) / sizeof(Object),
	};

	srand(12345);

	FILE* image_file = fopen("output.ppm", "w");
	fprintf(image_file, "P3\n%d %d\n255\n", image_width, image_height);
	for (int pixel_y = 0; pixel_y < image_height; pixel_y++) {
		/*printf("Rendering %d/%d", pixel_y + 1, image_height);*/
		for (int pixel_x = 0; pixel_x < image_width; pixel_x++) {
			vec3 total_color = { 0.0, 0.0, 0.0 };
			for (int sample = 0; sample < samples_per_pixel; sample++) {
				double offset_x = 1.0 * (frand() - 0.5);
				double offset_y = 1.0 * (frand() - 0.5);
				double fac_x = (pixel_x + 0.5 + offset_x) / (image_width - 1);
				double fac_y = (pixel_y + 0.5 + offset_y) / (image_height - 1);
				vec3 pixel_pos = {
					.x = (fac_x - 0.5) * viewport_width,
					.y = 0.0,
					.z = -(fac_y - 0.5) * viewport_height,
				};
				Ray ray = ray_between(eye_pos, pixel_pos);
				vec3 ray_color = trace_ray(ray, world, max_bounces);
				total_color = vec3_add(total_color, ray_color);
			}
			vec3 pixel_color = vec3_scale(total_color, 1.0 / samples_per_pixel);
			int red_byte = 255.999 * sqrt(pixel_color.x);
			int green_byte = 255.999 * sqrt(pixel_color.y);
			int blue_byte = 255.999 * sqrt(pixel_color.z);
			fprintf(image_file, "%d %d %d\n", red_byte, green_byte, blue_byte);
		}
		/*putchar('\n');*/
	}
	fclose(image_file);
	return 0;
}
