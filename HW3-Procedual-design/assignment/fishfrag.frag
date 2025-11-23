#version 300 es
precision mediump float;

// --- 标准 GLSL Canvas / ShaderToy 输入 ---
out vec4 colour_out;
uniform vec2 u_resolution; // 画布宽高
uniform float u_time;      // 时间 (原 vu_time)
uniform vec2 u_mouse;      // 鼠标位置 (保留接口，方便扩展)

// SDF Functions (Logic remains the same)
float sdEgg( in vec2 p, in float he, in float ra, in float rb, in float bu )
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

float getRuffleNoise(float angle, float time) {
    float noise = sin(angle * 8.0 + time * 1.5) * 0.5;
    noise += sin(angle * 17.0 - time * 2.0) * 0.25;
    noise += sin(angle * 37.0 + time * 3.5) * 0.125;
    return noise * 0.1; 
}

float sdBettaTailShape( vec2 p, vec2 c, float r, float iTime )
{
    float angle = atan(p.y, p.x);
    float ruffleOffset = getRuffleNoise(angle, iTime) * r;
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

float smoothUnion(float d1, float d2, float k) {
    float h = clamp(0.5 + 0.5 * (d2 - d1) / k, 0.0, 1.0);
    return mix(d2, d1, h) - k * h * (1.0 - h);
}

void main() {
    // --- 坐标系统转换 ---
    // 1. 获取归一化坐标 [0, 1] (左下角为0,0)
    vec2 st = gl_FragCoord.xy / u_resolution.xy;
    
    // 2. 线性变换到 [-1, 1] (中心为0,0)，以匹配你原本的 vPosition 逻辑
    // 公式: new_coord = old_coord * 2.0 - 1.0
    vec2 centcoord = st * 2.0 - 1.0;

    // (注意：这里没有进行纵横比修正，因为你的要求是保持你原本的数值逻辑
    // 如果画面被拉伸，你可能需要 p.x *= u_resolution.x / u_resolution.y)

    // --- 1. THE HEAD (ID:0) ---
    // Original Center: (0.650, 0.460) -> New: (0.65*2 - 1, 0.46*2 - 1) = (0.30, -0.08)
    vec2 center0 = centcoord - vec2(0.300, -0.080);
    
    // Dimensions are scaled by 2.0
    // a0 (0.32, 0.38) -> (0.64, 0.76)
    vec2 a0 = vec2(0.640, 0.760); 
    // b0 (0.07, 0.06) -> (0.14, 0.12)
    vec2 b0 = vec2(0.140, 0.120); 
    // w0 0.088 -> 0.176
    float w0 = 0.176;             
    float d0 = sdOrientedVesica(center0, a0, b0, w0);

    // --- 2. THE BODY (ID:1) ---
    // Original Center: (0.820, 0.670) -> New: (0.64, 0.34)
    vec2 center1 = centcoord - vec2(0.640, 0.340);
    
    // Rotation angles do NOT change
    float theta1 = -4.256;
    mat2 rotate1 = mat2(cos(theta1), -sin(theta1), sin(theta1), cos(theta1));
    
    // Dimensions scaled by 2.0
    float he1 = 0.488; // 0.244 * 2
    float ra1 = 0.204; // 0.102 * 2
    float rb1 = 0.144; // 0.072 * 2
    float bu1 = 0.528; // 0.264 * 2
    float d1 = sdEgg(rotate1 * center1, he1, ra1, rb1, bu1);

    // --- 3. THE TAIL (ID:2) ---
    // Original Center: (0.570, 0.570) -> New: (0.14, 0.14)
    vec2 center2 = centcoord - vec2(0.140, 0.140);
    
    float theta2 = -4.360;
    mat2 rotate2 = mat2(cos(theta2), -sin(theta2), sin(theta2), cos(theta2));
    
    // c2 is a direction vector (sin/cos), so it stays normalized (NO SCALE)
    vec2 c2 = vec2(0.920, 0.280); 
    
    // Radius scaled by 2.0
    float r2 = 1.000; // 0.500 * 2
    // 如果想要动起来，可以将 0.200 改为 u_time
    float itime2 = 0.200; 
    float d2 = sdBettaTailShape(rotate2 * center2, c2, r2, itime2);

    // --- 4. COMBINE ---
    // Smooth factor k also needs to be scaled by 2.0 to maintain relative softness
    // k = 0.048 -> 0.096
    float smooth_k = 0.096;

    float d_headbody = smoothUnion(d0, d1, smooth_k);
    float d_tail = d2;
    float d_temp = smoothUnion(d_headbody, d_tail, smooth_k);

    // Optional: Use smoothstep for anti-aliasing instead of hard step
    // 2.0 / 700.0 represents roughly 2 pixels of blur for smoothness
    float aa = 2.0 / 700.0; 
    float mask = 1.0 - smoothstep(0.0 - aa, 0.0 + aa, d_temp);
    
    // Output color
    // 这里使用 centcoord 代替了原来的 vPosition
    vec3 colorwoutAlpha = vec3(centcoord * 0.5 + 0.5, 1.0);
    colour_out = mask * vec4(colorwoutAlpha, 1.0);
}