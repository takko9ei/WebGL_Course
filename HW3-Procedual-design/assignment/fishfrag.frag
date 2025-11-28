#version 300 es
#define PI 3.14159265
#define E 2.71828

precision mediump float;

out vec4 colour_out;
uniform vec2 u_resolution;
uniform float u_time;
uniform vec2 u_mouse;

// Function Definitions

// SDF Functions
// Chinese comment use UTF-8

float remap(float value, float origFrom, float origTo, float targetFrom, float targetTo) {
    float normalizedValue = (value - origFrom) / (origTo - origFrom); // Normalize to [0, 1]
    return targetFrom + (targetTo - targetFrom) * normalizedValue; // Map to target range
}

float getRuffleNoise(float angle, float time,float noiseScale) {
    // numbers here has no special meaning, just random values
    float noise = sin(angle * 8.3 + time * 2.0) * 0.5;
    noise += sin(angle * 17.7 - time * 3.5) * 0.25;
    noise += sin(angle * 37.9 + time * 4.5) * 0.125;
    return noise * noiseScale; //the noise strength
}

vec3 getBodyPattern(vec2 uv, vec3 baseCol, vec3 scaleCol, vec3 hintCol, float scale, float rot, float iTime){
    // geo calc
    mat2 rotMat = mat2(cos(rot), -sin(rot), sin(rot), cos(rot));
    vec2 st = rotMat * uv;
    
    st *= scale * 7.0; 
    st.x += step(1.0, mod(st.y, 2.0)) * 0.5;
    
    vec2 cell_uv = fract(st);
    cell_uv.y = remap(cell_uv.y, 0.0, 1.0, 0.0, 0.65); // 将 y 方向压缩为等边三角形高度比例

    // dist
    float d = length(cell_uv - vec2(0.5, 0.0));
    
    float scale_shape = smoothstep(0.65, 0.45, d); 
    
    // thickness
    float shadow = smoothstep(0.2, 0.45, d);

    // color mix
    vec3 base_color = baseCol; 
    vec3 scale_color = scaleCol;
    
    vec3 iridescent = hintCol * sin(uv.y * 5.0 + iTime*5.0 * 0.5);

    // 混合：底色 + 鳞片形状 * (鳞片色 + 偏光)
    // 减去 shadow 让边缘变暗，而不是变亮，这样看起来像凸起而不是孔洞
    vec3 final_color = base_color + scale_shape * (scale_color + iridescent) - shadow * 0.2;

    return clamp(final_color, 0.0, 1.0);
}

vec3 getFinPattern(vec2 uv, vec3 mainCol, vec3 subCol, float density, float rot, vec2 st, float scale, float iTime){
    vec2 u0 = uv;
    float speed = 0.5;
    uv*=scale;
    mat2 rotMat = mat2(cos(rot), sin(rot), -sin(rot), cos(rot));
    uv = rotMat * uv;
    uv -= st;
    float shapeFactor = cos(clamp(uv.x * 1.2, -1.57, 1.57));
    uv.y /= (shapeFactor * shapeFactor + 0.1);
    for (int i = 0; i < 4; i++) {
              float fi = float(i);
              float offset = 0.15 / (fi + 1.0) * sin(fi * 3.0 * uv.y + iTime * speed + 0.3 * fi) + 0.5;
              uv.x += offset;
              uv.y += 0.1 / (fi + 1.0) * sin(fi * 4.0 * uv.x + iTime * speed * 0.8);
    }
    
    float pattern = sin(uv.x * density + uv.y * density);
    pattern = abs(pattern);
    pattern = 0.50 / max(pattern, 0.001); // max防止除0
    pattern = pow(pattern,1.8);
    pattern = min(pattern, 2.0);
    vec3 col = mainCol;
    col += vec3(0.2, 0.1, 0.4) * sin(u0.y * 2.0 + uv.y);
    col *= pattern;
    col += subCol * (1.0 - pattern);
    float mask = smoothstep(0.0, 0.8, shapeFactor);
    mask = pow(mask, 3.0);
    //float vFade = 1.0 - smoothstep(0.5, 1.0, abs(u0.y));

    col *= mask;// * vFade;
    return col;
}

float sdEgg( vec2 p, float he, float ra, float rb, float bu )
{
    float r = 0.5*(he + ra+rb)/bu;
    float da = r - ra;
    float db = r - rb;
    float y = (db*db - da*da - he*he)/(2.0*he);
    float x = sqrt(da*da - y*y);
    p.x = abs(p.x);
    float k = p.y*x - p.x*y;
    if( k>0.0 && k<he*(p.x+x) )
        return length(p+vec2(x,y))-r;
    return min( length(p)-ra,
                length(vec2(p.x,p.y-he))-rb );
}

float sdBettaTailShape( vec2 p, vec2 c, float r, float iTime )
{
    float angle = atan(p.y, p.x);
    float ruffleOffset = getRuffleNoise(angle, iTime,0.12) * r;
    p.x = abs(p.x);
    float l = length(p) - (r + ruffleOffset);
    float m = length(p-c*clamp(dot(p,c),0.0,r)); 
    return max(l,m*sign(c.y*p.x-c.x*p.y));
}

float sdOrientedVesica( vec2 p, vec2 a, vec2 b, float w )
{
    float r = 0.5*length(b-a);
    float d = 0.5*(r*r-w*w)/w;
    vec2 v = (b-a)/r;
    vec2 c = (b+a)*0.5;
    vec2 q = 0.5*abs(mat2(v.y,v.x,-v.x,v.y)*(p-c));
    vec3 h = (r*q.x<d*(q.y-r)) ? vec3(0.0,r,0.0) : vec3(-d,0.0,d+w);
    return length( q-h.xy) - h.z;
}

float sdMoon(vec2 p, float d, float ra, float rb, float intensity, float iTime )
{
    // --- 步骤 A: 预处理对称性 ---
    // 原函数利用了上下对称性。这意味着我们的噪声也会随之上下对称（镜像）。
    // 如果想要非对称的噪声，需要修改这行代码的逻辑，但为了保持 SDF 算法的稳定性，
    // 我们暂时保留这个对称性，这意味着月亮咬痕的上下两半是镜像的。
    p.y = abs(p.y);

    // --- 步骤 B: 计算相对于“小圆圆心”的角度 ---
    // 小圆圆心在 (d, 0)。
    // 我们需要计算当前点 p 到 (d, 0) 的向量。
    vec2 relP = p - vec2(d, 0.0);
    
    // 获取角度。atan(y, x) 返回弧度。
    float angle = atan(relP.y, relP.x);

    // --- 步骤 C: 计算扰动后的半径 ---
    
    // 这里的 rb 变成了变量。
    // 注意：我们把 getNoise 的结果加到 rb 上
    float noisy_rb = rb + getRuffleNoise(angle,iTime, intensity);

    // --- 步骤 D: 代入原公式 ---
    // 关键点：原本公式中所有用到 rb 的地方，现在全部都要换成 noisy_rb
    // 这样才能保证月亮的尖角 (intersection points) 计算是基于变形后的圆的
    
    float a = (ra*ra - noisy_rb*noisy_rb + d*d)/(2.0*d);
    float b = sqrt(max(ra*ra-a*a,0.0));
    
    if( d*(p.x*b-p.y*a) > d*d*max(b-p.y,0.0) )
          return length(p-vec2(a,b));
    
    return max( (length(p)-ra),
               -(length(p-vec2(d,0))-noisy_rb));
}

float sMinCubic( float a, float b, float k )
{
    float h = max( k - abs(a-b), 0.0 )/k;
    return min( a, b ) - h*h*h*k*(1.0/6.0);
}

float sDiffCubic( float d1, float d2, float k )
{
    float h = max( k - abs(d1 + d2), 0.0 )/k;
    return max( d1, -d2 ) + h*h*h*k*(1.0/6.0);
}

//input: 1d distance and 3d color. 
//this function returns the smooth union of d(sdf result) and color
vec4 sUnion(vec4 d1_c1, vec4 d2_c2, float k) {
    float h = clamp( 0.5 + 0.5 * (d1_c1.x - d2_c2.x) / k, 0.0, 1.0 );
    float dist = sMinCubic(d1_c1.x, d2_c2.x, k);
    vec3 col = mix(d1_c1.yzw, d2_c2.yzw, h);
    return vec4(dist, col);
}

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

vec3 getWavePattern(vec2 uv, vec3 mainCol, vec3 subCol, float scale, float iTime) {
    vec2 scaledUV = uv * scale * 6.0; // Scale to control getWavePattern size

    // Get distances to the two nearest cell centers
    vec2 dists = voronoi(scaledUV, iTime);

    // The first distance gives a nice gradient within the cell
    float cell_interior = dists.x;
    
    // The difference between the two distances defines the cell edge
    float edge_dist = dists.y - dists.x;

    // --- Coloring ---
    // Base color for the water, darker blue
    vec3 water_color = mainCol;
    // Make the center of the getWavePattern polygons a bit brighter
    water_color = mix(subCol, water_color, smoothstep(0.0, 0.4, cell_interior));
    
    // --- Sparkle & Bloom (with AA and Smoother Shimmer) ---

    // 1. Smoother, organic flicker mask
    float flicker_noise_sample = noise(scaledUV * 0.5 + iTime * 0.2);
    // Use a second, smoother noise for shimmer instead of a raw hash
    float shimmer_noise = noise(scaledUV * 1.5 + iTime * -0.5);
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
// --- New Function: Water Ripple Effect based on your reference ---
// Adapt the HLSL logic to GLSL
vec2 getWaterRipple(vec2 uv, vec2 center, float time) {
    // 1. Calculate vector from current pixel to the ripple center
    // HLSL: float2 dv = center - i.uv;
    vec2 dv = center - uv;

    // 2. Calculate distance
    // HLSL: float dis = sqrt(...)
    float dis = length(dv);

    // 3. Ripple parameters (You can tweak these)
    float frequency = 15.0;  // Corresponds to _Ctor (How dense the rings are)
    float speed = 5.0;       // Corresponds to _timeCtor (How fast they move)
    float maxRadius = 1.2;   // Corresponds to _max_dis (Ripple range)
    float amplitude = 0.03;  // Strength of the distortion

    // 4. Calculate Sine Factor (The wave)
    // HLSL: sin(dis * _Ctor + _Time.y * _timeCtor)
    // Note: Usually "- time" looks like expanding waves, "+ time" looks like imploding.
    float sinFactor = sin(dis * frequency - time * speed);

    // 5. Calculate Damping/Falloff
    // HLSL: max(0, _max_dis - dis)
    // Limits the effect to a specific radius and fades it out at the edges
    float damp = max(0.0, maxRadius - dis);
    
    // Optional: Add extra smoothness to the edge
    damp = pow(damp, 2.0); 

    // 6. Calculate Direction
    // HLSL: normalize(dv)
    vec2 dir = normalize(dv);
    
    // Safety check to avoid NaN when uv is exactly at center
    if (dis < 0.0001) dir = vec2(0.0); 

    // 7. Final Offset
    // HLSL: dv1 * sinFactor * max(...) * step(...)
    // We multiply by 'amplitude' to control strength
    return dir * sinFactor * damp * amplitude;
}

void main() {
    vec2 st = gl_FragCoord.xy / u_resolution.xy;

    // st: left bottom (0,0) to right top (1,1), linear
    // we need x be distorted, y be linear
    //用x的x坐标控制x的拉伸
    float xScaleFactor = remap(sin(u_time),-1.0,1.0,0.0,2.0);
    //st.x*=xScaleFactor;
    
    vec2 centcoord = st * 2.0 - 1.0;

    vec2 mouse_uv = u_mouse.xy / u_resolution.xy;
    vec2 rippleCenter = mouse_uv * 2.0 - 1.0;
    
    // Option B: Fixed center (e.g., (0,0)) if mouse is not active
    if (u_mouse.x < 10.0) { rippleCenter = vec2(0.0, 0.0); }

    // Calculate the offset
    vec2 rippleOffset = getWaterRipple(centcoord, rippleCenter, u_time);

    // Apply the offset to the coordinate system
    centcoord += rippleOffset;

    // centcoord: left bottom (-1,-1) to right top (1,1), linear

    vec3 color_body = getBodyPattern(centcoord,vec3(0.4314, 0.1765, 0.502), vec3(0.2, 0.7, 0.8), vec3(0.3765, 0.1608, 0.8039),1.952,2.032, u_time); 
    vec3 color_backfin = getFinPattern(centcoord,vec3(0.1608, 0.5961, 0.7843),vec3(0.3451, 0.0235, 0.0235),23.0,-0.192,vec2(-0.200,0.590),1.128,u_time);
    vec3 color_bellyfin = getFinPattern(centcoord,vec3(0.1647, 0.5686, 0.8784),vec3(0.3451, 0.0235, 0.0235),21.0,-0.856,vec2(0.480,0.630),1.80,u_time);
    vec3 color_tail = getFinPattern(centcoord,vec3(0.1412, 0.502, 0.7804),vec3(0.3451, 0.0235, 0.0235),18.0,-0.416,vec2(-0.790,0.580),1.0,u_time);
    vec3 color_water = getWavePattern(centcoord, vec3(0.3686, 0.7176, 0.7686), vec3(0.549, 0.7922, 0.8039), 0.7, u_time);
    vec2 uv = centcoord;
    // --- 1. THE HEAD (ID:0) ---
    vec2 p0 = centcoord - vec2(0.300, -0.080);
    vec2 a0 = vec2(0.640, 0.760); 
    vec2 b0 = vec2(0.140, 0.120); 
    float w0 = 0.176;             
    float d0 = sdOrientedVesica(p0, a0, b0, w0);
    // --- 2. THE BODY (ID:1) ---
    vec2 p1 = centcoord - vec2(0.640, 0.340);
    float theta1 = -4.256;
    mat2 rotate1 = mat2(cos(theta1), -sin(theta1), sin(theta1), cos(theta1));
    float he1 = 0.488;
    float ra1 = 0.204;
    float rb1 = 0.144;
    float bu1 = 0.528;
    float d1 = sdEgg(rotate1 * p1, he1, ra1, rb1, bu1);
    // --- 3. THE TAIL (ID:2) ---
    vec2 p2 = centcoord - vec2(0.140, 0.140);
    float theta2 = -4.360;
    mat2 rotate2 = mat2(cos(theta2), -sin(theta2), sin(theta2), cos(theta2));
    vec2 c2 = vec2(0.920, 0.280); 
    float r2 = 1.00; 
    float itime2 = u_time; 
    float d2 = sdBettaTailShape(rotate2 * p2, c2, r2, itime2);
   //--- 4. THE BACK FIN (ID:3) ---
   vec2 p3 = centcoord - vec2(-0.610,0.410);
   float theta3 = 0.448;
   mat2 rotate3 = mat2(cos(theta3), -sin(theta3), sin(theta3), cos(theta3));
   float cfactor = 0.736;
   mat2 compress = mat2(cfactor*st.x,0,0,1);
   float pm3 = -0.440;
   float ra3 = 0.480;
   float rb3 = 0.688;
   float iTime3 = u_time;
   float intensity3 = 0.042;
   float d3 = sdMoon(rotate3*compress*p3,pm3,ra3,rb3,intensity3,iTime3);
	// --- 5. The belly fin(ID:4) ---
   vec2 p4 = centcoord - vec2(0.390,0.050);
   float theta4 = -5.408;
   mat2 rotate4 = mat2(cos(theta4), -sin(theta4), sin(theta4), cos(theta4));
   float pm4 = -0.752;
   float ra4 = 0.376;
   float rb4 = 0.616;
   float iTime4 = u_time;
   float intensity4 = 0.042;
   float d4 = sdMoon(rotate4*p4,pm4,ra4,rb4,intensity4,iTime4);    

    // COLOR

    // --- 6. CLIP AND COMBINE ---
    float smooth_k = 0.096;
    float clip_margin = -0.01;
    float fin_k = 0.04;
    
    float d_headbody = sMinCubic(d0, d1, smooth_k);
    float d_tail = d2;
    float d_backfin = d3;
    float d_bellyfin = d4;

    float d_body_core = d_headbody + clip_margin;
    float d_tail_core = d_tail + clip_margin-0.1;

    d_backfin = sDiffCubic(d_backfin, d_body_core, 0.02);
    d_bellyfin = sDiffCubic(d_bellyfin, d_body_core, 0.02);
    d_bellyfin = sDiffCubic(d_bellyfin, d3, 0.02);
    d_bellyfin = sDiffCubic(d_bellyfin,d_tail_core,0.02);

    vec4 current_res = vec4(d_headbody, color_body);
    current_res = sUnion(current_res, vec4(d_backfin, color_backfin), fin_k);
    current_res = sUnion(current_res, vec4(d_bellyfin, color_bellyfin), fin_k);

    if(d_backfin<0.0&&d_tail<0.01){
        current_res.yzw = color_tail;
    }

    current_res = sUnion(current_res,vec4(d_tail,color_tail), fin_k);

    float d_final = current_res.x;

    float mask = step(d_final,0.0);

    float mask_eye = step(length(centcoord - vec2(0.800, 0.53)) , 0.0450);
    mask_eye = remap(mask_eye, 0.0, 1.0, 1.0, 0.0);
    float mask_eye_highlight = step(length(centcoord - vec2(0.790, 0.54)) , 0.01);
    //mask_eye_highlight = remap(mask_eye_highlight, 0.0, 1.0, 1.0, 0.0);
    //temp cancel
    
    // Output color
    //colour_out = vec4(current_res.yzw,1.0);
    colour_out = mask * vec4(current_res.yzw, 1.0)+ (1.0 - mask) * vec4(color_water, 1.0);
    //the most violent way to draw eye
    colour_out.xyz = colour_out.xyz * mask_eye;
    colour_out.xyz += mask_eye_highlight * vec3(1.0, 1.0, 1.0);
    colour_out = clamp(colour_out, 0.0, 1.0);
}