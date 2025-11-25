#version 300 es
precision mediump float;

out vec4 colour_out;

uniform vec2 u_resolution;
uniform float u_time;
uniform vec2 u_mouse;

// Function to generate a fish scale pattern value for a given UV coordinate.
// The geometry is static, but the returned value can be used for dynamic coloring.
float scale_pattern(vec2 uv) {
    // Scale UV to create more scales
    uv *= 15.0;

    // Offset every other row to create the overlapping pattern
    uv.x += step(1.0, mod(uv.y, 2.0)) * 0.5;

    // Get the fractional part of the UV to work within a single grid cell
    vec2 cell_uv = fract(uv);

    // Calculate distance from the bottom-center of the cell
    // This creates an arc, which looks like a scale
    float dist_from_center = length(cell_uv - vec2(0.5, 0.0));

    return dist_from_center;
}

void main() {
    // Normalize fragment coordinates and correct aspect ratio
    vec2 st = gl_FragCoord.xy / u_resolution.xy;
    st.x *= u_resolution.x / u_resolution.y;

    // --- 1. Geometry ---
    // Get the basic geometric value for the scales. This is static.
    float d = scale_pattern(st);

    // --- 2. Color and Glow ---
    
    // Create the glowing edge effect, similar to getFinPattern.
    // We use sin on the distance value. Where the distance creates a nice gradient (like 0, 1, 2...),
    // sin will create waves. The 1/abs(sin) trick creates bright lines at the zero-crossings.
    float glow_pattern = 0.04 / max(abs(sin(d * 6.28)), 0.001); // 6.28 is 2*PI
    glow_pattern = clamp(glow_pattern, 0.0, 1.0);

    // Define a base color for the scales
    vec3 base_color = vec3(0.1, 0.5, 0.7); // A bluish-green

    // Create an iridescent, shimmering effect that moves over time.
    // This is inspired by: col += vec3(0.2, 0.1, 0.4) * sin(u0.y * 2.0 + uv.y);
    // We add u_time to make the color bands move vertically.
    vec3 iridescent_color = vec3(0.3, 0.2, 0.5) * sin(st.y * 8.0 + u_time * 0.5);

    // Combine base color and iridescent shimmer
    vec3 final_color = base_color + iridescent_color;

    // Apply the glowing edge pattern, same as `col *= pattern;` in getFinPattern
    final_color *= glow_pattern;

    // --- 4. Final Output ---
    colour_out = vec4(final_color, 1.0);
}
