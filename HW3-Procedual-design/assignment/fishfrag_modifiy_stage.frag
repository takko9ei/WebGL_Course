#version 300 es
precision mediump float;

out vec4 colour_out;
uniform vec2 u_resolution;
uniform float u_time;
uniform vec2 u_mouse;

// SDF Functions
// Chinese comment use UTF-8
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

float getRuffleNoise(float angle, float time,float noiseScale) {
    float noise = sin(angle * 8.3 + time * 2.0) * 0.5;
    noise += sin(angle * 17.7 - time * 3.5) * 0.25;
    noise += sin(angle * 37.9 + time * 4.5) * 0.125;
    return noise * noiseScale; 
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
    p.y = abs(p.y);
    vec2 relP = p - vec2(d, 0.0);
    float angle = atan(relP.y, relP.x);
    float noisy_rb = rb + getRuffleNoise(angle,iTime, intensity);
    float a = (ra*ra - noisy_rb*noisy_rb + d*d)/(2.0*d);
    float b = sqrt(max(ra*ra-a*a,0.0));
    
    // 如果 d 为负数，这里的逻辑可能会出错，导致内部返回正距离
    if( d*(p.x*b-p.y*a) > d*d*max(b-p.y,0.0) )
          return length(p-vec2(a,b));
    return max( (length(p)-ra),
               -(length(p-vec2(d,0))-noisy_rb));
}

float smin_cubic( float a, float b, float k )
{
    float h = max( k - abs(a-b), 0.0 )/k;
    return min( a, b ) - h*h*h*k*(1.0/6.0);
}

// 混合函数：同时混合距离(x)和颜色(yzw)
vec4 opSmoothUnion(vec4 d1_c1, vec4 d2_c2, float k) {
    float h = clamp( 0.5 + 0.5 * (d1_c1.x - d2_c2.x) / k, 0.0, 1.0 );
    float dist = smin_cubic(d1_c1.x, d2_c2.x, k);
    vec3 col = mix(d1_c1.yzw, d2_c2.yzw, h);
    return vec4(dist, col);
}

void main() {
    vec2 st = gl_FragCoord.xy / u_resolution.xy;
    vec2 centcoord = st * 2.0 - 1.0;

    // --- Colors ---
    vec3 color_body     = vec3(1.0, 0.4, 0.2); 
    vec3 color_backfin  = vec3(0.8, 0.1, 0.1); 
    vec3 color_bellyfin = vec3(0.2, 0.6, 1.0); 
    vec3 color_tail     = vec3(0.7, 0.2, 0.6); 

    // --- 1. HEAD ---
    vec2 center0 = centcoord - vec2(0.300, -0.080);
    vec2 a0 = vec2(0.640, 0.760); 
    vec2 b0 = vec2(0.140, 0.120); 
    float w0 = 0.176;             
    float d0 = sdOrientedVesica(center0, a0, b0, w0);

    // --- 2. BODY ---
    vec2 center1 = centcoord - vec2(0.640, 0.340);
    float theta1 = -4.256;
    mat2 rotate1 = mat2(cos(theta1), -sin(theta1), sin(theta1), cos(theta1));
    float he1 = 0.488;
    float ra1 = 0.204;
    float rb1 = 0.144;
    float bu1 = 0.528;
    float d1 = sdEgg(rotate1 * center1, he1, ra1, rb1, bu1);

    // --- 3. TAIL ---
    vec2 center2 = centcoord - vec2(0.140, 0.140);
    float theta2 = -4.360;
    mat2 rotate2 = mat2(cos(theta2), -sin(theta2), sin(theta2), cos(theta2));
    vec2 c2 = vec2(0.920, 0.280); 
    float r2 = 1.00; 
    float d2 = sdBettaTailShape(rotate2 * center2, c2, r2, u_time);

    // --- 4. BACK FIN ---
    vec2 center3 = centcoord - vec2(-0.610,0.410);
    float theta3 = 0.448;
    mat2 rotate3 = mat2(cos(theta3), -sin(theta3), sin(theta3), cos(theta3));
    float cfactor = 0.736;
    mat2 compress = mat2(cfactor*st.x,0,0,1);
    float pm3 = -0.440; 
    float ra3 = 0.480;
    float rb3 = 0.688;
    float intensity3 = 0.042;
    float d3 = sdMoon(rotate3*compress*center3,pm3,ra3,rb3,intensity3,u_time);

    // --- 5. BELLY FIN ---
    vec2 center4 = centcoord - vec2(0.390,0.050);
    float theta4 = -5.408;
    mat2 rotate4 = mat2(cos(theta4), -sin(theta4), sin(theta4), cos(theta4));
    float pm4 = -0.752;
    float ra4 = 0.376;
    float rb4 = 0.616;
    float intensity4 = 0.042;
    float d4 = sdMoon(rotate4*center4,pm4,ra4,rb4,intensity4,u_time);

    // --- 6. COMBINE ---
    float smooth_k = 0.096;
    float fin_k = 0.04;

    // 1. Head (d0) + Body (d1) -> BodyGroup
    // 头部和身体颜色相同，所以不需要计算复杂的颜色插值，直接取 min 即可，但为了形状平滑沿用 smin
    float d_headbody = smin_cubic(d0, d1, smooth_k);
    vec3 col_current = color_body; 

    // 2. + Backfin (d3)
    // 计算颜色混合权重 h_back: 当 d3 (背鳍) 越小于 d_headbody (更内部) 时，h_back 越趋近 1.0
    float h_back = clamp(0.5 + 0.5 * (d_headbody - d3) / fin_k, 0.0, 1.0);
    float d_headbodyfin = smin_cubic(d_headbody, d3, fin_k);
    // 颜色插值
    col_current = mix(col_current, color_backfin, h_back);

    // 3. + Bellyfin (d4)
    float h_belly = clamp(0.5 + 0.5 * (d_headbodyfin - d4) / fin_k, 0.0, 1.0);
    d_headbodyfin = smin_cubic(d_headbodyfin, d4, fin_k);
    col_current = mix(col_current, color_bellyfin, h_belly);

    // 4. + Tail (d2)
    // [关键逻辑] 动态 K 值
    // 我们依据 d3 (背鳍距离) 来调整尾巴的融合平滑度。
    // smoothstep(-0.02, 0.1, d3): 
    //   当 d3 < -0.02 (在背鳍内部) -> 结果为 0.0 -> k_tail 变为 0.0 (硬边)
    //   当 d3 > 0.1 (远离背鳍，即在身体附近) -> 结果为 1.0 -> k_tail 变为 smooth_k (平滑)
    float k_tail_factor = smoothstep(-0.02, 0.1, d3); 
    float dynamic_k = smooth_k * k_tail_factor;

    // 防止 dynamic_k 过小导致除零错误，做个小判断
    float d_final;
    float h_tail;
    
    if(dynamic_k < 0.001) {
        // 硬连接 (Min)
        d_final = min(d_headbodyfin, d2);
        // 硬连接时的颜色切换：谁小(谁在内部)显谁
        h_tail = (d_headbodyfin < d2) ? 0.0 : 1.0; 
    } else {
        // 平滑连接 (Smin)
        h_tail = clamp(0.5 + 0.5 * (d_headbodyfin - d2) / dynamic_k, 0.0, 1.0);
        d_final = smin_cubic(d_headbodyfin, d2, dynamic_k);
    }
    
    col_current = mix(col_current, color_tail, h_tail);

    // [关键逻辑] 强制相交处颜色
    // 规则：Backfin (d3) 与 Tail (d2) 相交处，必须是 Backfin 的颜色
    if (d3 < 0.0 && d2 < 0.0) {
        col_current = color_backfin;
    }

    float mask = step(d_final, 0.0);
    
    // Output color
    // 使用计算好的 col_current 替换原来的坐标颜色
    colour_out = mask * vec4(col_current, 1.0);
}