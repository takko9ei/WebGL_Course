#version 300 es
precision mediump float;

out vec4 colour_out;

uniform vec2 u_resolution;
uniform float u_time;

// Provides a pseudo-random value for a 2D coordinate
float hash( vec2 p ) {
    p = vec2( dot(p,vec2(127.1,311.7)),
              dot(p,vec2(269.5,183.3)) );
    return fract(sin(dot(p, vec2(12.9898, 78.233))) * 43758.5453);
}

// Smoothly interpolated noise function (Value Noise)
float noise( vec2 x ) {
    vec2 i = floor(x);
    vec2 f = fract(x);

    // Smoothstep for interpolation
    vec2 u = f*f*(3.0-2.0*f);

    // Get random values for the 4 corners of the cell
    float a = hash( i + vec2(0.0,0.0) );
    float b = hash( i + vec2(1.0,0.0) );
    float c = hash( i + vec2(0.0,1.0) );
    float d = hash( i + vec2(1.0,1.0) );

    // Bilinearly interpolate
    return mix(mix( a, b, u.x), mix( c, d, u.x), u.y);
}

// Returns the distances to the nearest and second-nearest cell centers
vec2 voronoi(vec2 uv, float time) {
    vec2 g = floor(uv); // Grid cell
    vec2 f = fract(uv); // Position within cell
    
    vec2 closest = vec2(1.0); // .x is dist to nearest, .y is dist to 2nd nearest

    // Check 3x3 grid of cells around the current one
    for (int y = -1; y <= 1; y++) {
        for (int x = -1; x <= 1; x++) {
            vec2 offset = vec2(x, y);
            vec2 cell_pos = g + offset;
            
            // Animate each cell point pseudo-randomly
            float h = hash(cell_pos);
            vec2 point_pos = 0.5 + 0.4 * vec2(sin(time * 0.7 + h * 6.28), cos(time * 0.7 + h * 6.28)); // 2*PI
            
            float dist = length(offset + point_pos - f);

            // Update the two closest distances
            if (dist < closest.x) {
                closest.y = closest.x;
                closest.x = dist;
            } else if (dist < closest.y) {
                closest.y = dist;
            }
        }
    }
    return closest;
}

vec3 wave(vec2 uv, float scale, float iTime) {
    vec2 scaled_uv = uv * scale * 6.0; // Scale to control wave size

    // Get distances to the two nearest cell centers
    vec2 dists = voronoi(scaled_uv, iTime);

    // The first distance gives a nice gradient within the cell
    float cell_interior = dists.x;
    
    // The difference between the two distances defines the cell edge
    float edge_dist = dists.y - dists.x;

    // --- Coloring ---
    // Base color for the water, darker blue
    vec3 water_color = vec3(0.0, 0.2, 0.5);
    // Make the center of the wave polygons a bit brighter
    water_color = mix(vec3(0.1, 0.4, 0.7), water_color, smoothstep(0.0, 0.4, cell_interior));
    
    // --- Sparkle & Bloom (with AA and Smoother Shimmer) ---

    // 1. Smoother, organic flicker mask
    float flicker_noise_sample = noise(scaled_uv * 0.5 + iTime * 0.2);
    // Use a second, smoother noise for shimmer instead of a raw hash
    float shimmer_noise = noise(scaled_uv * 1.5 + iTime * -0.5);
    float shimmer = (0.6 + 0.4 * sin(iTime * 8.0 + shimmer_noise * 6.28));
    
    float flicker_mask = mix(0.5, 1.0, flicker_noise_sample);
    flicker_mask *= shimmer;

    // 2. Anti-aliased Bloom effect
    float edge_width = fwidth(edge_dist) * 1.5; // Get pixel width for AA
    // Sharp core with AA
    float sharp_glow = pow(1.0 - smoothstep(0.0, edge_width * 2.0, edge_dist), 40.0);
    // Softer bloom with AA
    float soft_glow = pow(1.0 - smoothstep(0.0, edge_width * 8.0, edge_dist), 10.0);

    // Combine the bloom layers
    float total_bloom = sharp_glow * 1.0 + soft_glow * 0.3;

    // 3. Combine bloom with the flicker mask
    float sparkle = total_bloom * flicker_mask;
    
    vec3 sparkle_color = vec3(0.9, 1.0, 1.0);
    
    // Combine water color with the sparkling highlights
    vec3 final_color = water_color + sparkle * sparkle_color;

    return clamp(final_color, 0.0, 1.0);
}

void main() {
    // Normalize fragment coordinates. No aspect ratio correction as requested.
    vec2 uv = gl_FragCoord.xy / u_resolution.xy;

    // Define scale and time
    float scale = 1.0;
    float iTime = u_time;

    // Get the wave color
    vec3 final_color = wave(uv, scale, iTime);

    // Output final color
    colour_out = vec4(final_color, 1.0);
}
