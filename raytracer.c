#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

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

vec3 vec3_random(void) {
	return vec3_normalize((vec3) {
		.x = frand(),
		.y = frand(),
		.z = frand(),
	});
}

/*static const vec3 UP = { 0.0, 0.0, 1.0 };*/
/*static const vec3 RIGHT = { 0.0, 0.0, -1.0 };*/
/*static const vec3 FORWARD = { -1.0, 0.0, 0.0 };*/


/*static const vec3 RED = { 1.0, 0.0, 0.0 };*/
/*static const vec3 GREEN = { 0.0, 1.0, 0.0 };*/
/*static const vec3 BLUE = { 0.0, 0.0, 1.0 };*/
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
	bool outfacing;
	double t;
	vec3 point;
	vec3 normal;
} Hit;

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
		hit->outfacing = (vec3_dot(ray.direction, hit->normal) < 0.0);
		hit->normal = hit->outfacing ? hit->normal : vec3_scale(hit->normal, -1.0);
	}
	return true;
}


vec3 get_world_color(Ray ray) {
	/*const vec3 sun_dir = vec3_normalize((vec3) { -1.0, -1.0, -1.0 });*/
	/*float fac = fmax(0.0, vec3_dot(ray.dir, vec3_scale(sun_dir, -1.0)));*/
	static const vec3 sky_color = { 0.2, 0.4, 0.8 };
	double fac = 0.5 * (ray.direction.z + 1.0);
	return vec3_lerp(sky_color, WHITE, fac);
}

vec3 scatter_diffuse(vec3 incoming, vec3 normal) {
	vec3 outgoing = vec3_random();
	if (vec3_dot(outgoing, normal) <= 0.0) {
		outgoing = vec3_scale(outgoing, -1.0);
	}
	return outgoing;
}

vec3 lambertian_diffuse(vec3 incoming, vec3 normal, vec3 albedo) {
	const vec3 light_direction = vec3_normalize((vec3) { -1.0, -1.0, -1.0 });
	double intensity = fmax(0.0, vec3_dot(normal, vec3_scale(light_direction, -1.0)));
	vec3 color = (vec3) { 
		albedo.x * intensity, 
		albedo.y * intensity, 
		albedo.z * intensity 
	};
	return color;
}

bool ray_hit_world(Ray ray, double t_min, double t_max, Hit* hit) {

	static const Sphere sphere = {
		.center = { 0.0, 1.5, 0.0 },
		.radius = 0.5,
	};

	static const Sphere ground = {
		.center = { 0.0, 1.5, -100.5 },
		.radius = 100.0,
	};

	bool hit_anything = false;
	Hit hit_temp;
	if (ray_hit_sphere(ray, sphere, t_min, t_max, &hit_temp)) {
		hit_anything = true;
		t_max = hit_temp.t;
		if (hit != NULL) {
			*hit = hit_temp;
		}
	}
	if (ray_hit_sphere(ray, ground, t_min, t_max, &hit_temp)) {
		hit_anything = true;
		t_max = hit_temp.t;
		if (hit != NULL) {
			*hit = hit_temp;
		}
	}

        return hit_anything;
}

vec3 trace_ray(Ray ray, int bounces) {
	static const vec3 sphere_color = { 1.0, 1.0, 1.0 };
	Hit hit;
	if (bounces > 0 && ray_hit_world(ray, 0.0, INFINITY, &hit)) {
		vec3 diffuse_color = lambertian_diffuse(ray.direction, hit.normal, sphere_color);
		Ray bounced_ray = {
			.origin = hit.point,
			.direction = scatter_diffuse(ray.direction, hit.normal),
		};
		vec3 bounced_color = trace_ray(bounced_ray, bounces - 1);
		return vec3_lerp(diffuse_color, bounced_color, 0.5);
	}
	return get_world_color(ray);
}

int main(void) {

	static const int image_width = 512;
	static const int image_height = 512;
	static const double viewport_width = 1.0;
	static const double viewport_height = 1.0;
	static const double focal_length = 1.0;
	static const vec3 eye_pos = { 0.0, -focal_length, 0.0 };
	static const int max_bounces = 4;
	static const int samples_per_pixel = 8;

	srand(123456);

	FILE* image_file = fopen("output.ppm", "w");
	fprintf(image_file, "P3\n%d %d\n255\n", image_width, image_height);
	for (int pixel_y = 0; pixel_y < image_height; pixel_y++) {
		/*printf("Rendering %d/%d", pixel_y + 1, image_height);*/
		for (int pixel_x = 0; pixel_x < image_height; pixel_x++) {
			vec3 total_color = { 0.0, 0.0, 0.0 };
			for (int sample = 0; sample < samples_per_pixel; sample++) {
				double offset_x = frand() - 0.5;
				double offset_y = frand() - 0.5;
				double fac_x = ((double)pixel_x + 0.5 + offset_x) / (image_width - 1);
				double fac_y = ((double)pixel_y + 0.5 + offset_y) / (image_height - 1);
				vec3 pixel_pos = {
					.x = (fac_x - 0.5) * viewport_width,
					.y = 0.0,
					.z = -(fac_y - 0.5) * viewport_height,
				};
				Ray ray = ray_between(eye_pos, pixel_pos);
				vec3 sample_color = trace_ray(ray, max_bounces);
				total_color = vec3_add(total_color, sample_color);
			}
			vec3 pixel_color = vec3_scale(total_color, 1.0 / samples_per_pixel);
			int red_byte = 255.999 * pixel_color.x;
			int green_byte = 255.999 * pixel_color.y;
			int blue_byte = 255.999 * pixel_color.z;
			fprintf(image_file, "%d %d %d\n", red_byte, green_byte, blue_byte);
		}
		putchar('\n');
	}
	fclose(image_file);
	return 0;
}
