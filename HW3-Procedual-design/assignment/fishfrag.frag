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

float getRuffleNoise(float angle, float time) {
    float noise = sin(angle * 8.0 + time * 1.5) * 0.5;
    noise += sin(angle * 17.0 - time * 2.0) * 0.25;
    noise += sin(angle * 37.0 + time * 3.5) * 0.125;
    return noise * 0.12; //the noise strength
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

float getSideWobble(float dist, float time) {
    float noise = sin(dist * 10.0 - time * 3.0) * 0.5;
    noise += sin(dist * 20.0 + time * 4.0) * 0.30;
    
    // 关键点：乘以 dist * 0.15
    // 这样根部 (dist=0) 完全不动，越往末端摆动幅度越大，符合生物学
    return noise * 0.25 * dist; 
}

float sdBettaTailShape3edges(vec2 p, vec2 c, float r, float iTime) {
    float angle = atan(p.y, p.x);
    float ruffleOffset = getRuffleNoise(angle, iTime) * r;
    
    // 计算到圆弧的距离 'l'
    // 这里使用原始坐标 p 的长度，减去带噪声的半径
    float l = length(p) - (r + ruffleOffset);

    // 处理直边
    // 我们需要旋转坐标 p，来模拟侧边的摆动
    float dist = length(p); // 获取当前点离中心的距离
    float wobbleAngle = getSideWobble(dist, iTime); // 根据距离计算扭曲角度

    // 构建 2D 旋转矩阵
    float co = cos(wobbleAngle);
    float si = sin(wobbleAngle);
    mat2 rot = mat2(co, -si, si, co);

    // 旋转 p 得到专门用于计算侧边的临时坐标 pSide
    vec2 pSide = rot * p;

    // 标准扇形计算
    // 这里使用 pSide 来计算 m (侧边距离)，而不是原始 p
    pSide.x = abs(pSide.x);
    
    // 计算到直边的距离 'm'
    // 为了让侧边和圆弧完美衔接，clamp 的上限建议也加上 ruffleOffset (或者取近似值 r)
    // 但注意 ruffleOffset 是基于原始角度的
    float m = length(pSide - c * clamp(dot(pSide, c), 0.0, r));

    // 组合
    // 结合圆弧距离 l 和侧边距离 m
    // 这里的 sign 判断也必须使用旋转后的 pSide，否则内部填充会错位
    return max(l, m * sign(c.y * pSide.x - c.x * pSide.y));
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

float smin_cubic( float a, float b, float k )
{
    float h = max( k - abs(a-b), 0.0 )/k;
    return min( a, b ) - h*h*h*k*(1.0/6.0);
}

void main() {
    vec2 st = gl_FragCoord.xy / u_resolution.xy;
    
    vec2 centcoord = st * 2.0 - 1.0;

    // --- 1. THE HEAD (ID:0) ---
    vec2 center0 = centcoord - vec2(0.300, -0.080);
    vec2 a0 = vec2(0.640, 0.760); 
    vec2 b0 = vec2(0.140, 0.120); 
    float w0 = 0.176;             
    float d0 = sdOrientedVesica(center0, a0, b0, w0);

    // --- 2. THE BODY (ID:1) ---
    vec2 center1 = centcoord - vec2(0.640, 0.340);
    float theta1 = -4.256;
    mat2 rotate1 = mat2(cos(theta1), -sin(theta1), sin(theta1), cos(theta1));
    float he1 = 0.488;
    float ra1 = 0.204;
    float rb1 = 0.144;
    float bu1 = 0.528;
    float d1 = sdEgg(rotate1 * center1, he1, ra1, rb1, bu1);

    // --- 3. THE TAIL (ID:2) ---
    vec2 center2 = centcoord - vec2(0.140, 0.140);
    float theta2 = -4.360;
    mat2 rotate2 = mat2(cos(theta2), -sin(theta2), sin(theta2), cos(theta2));
    vec2 c2 = vec2(0.920, 0.280); 
    float r2 = 1.00; 
    float itime2 = u_time; 
    float d2 = sdBettaTailShape(rotate2 * center2, c2, r2, itime2);


    //--- 4. THE BACK FIN (ID:3) ---
    vec2 center30 = centcoord - vec2(0.210,0.310);
    float theta30 = 0.800;
    mat2 rotate30= mat2(cos(theta30), -sin(theta30), sin(theta30), cos(theta30));
    vec2 c30 = vec2(0.270,0.390); 
    float r30 = 0.5; 
    float itime30 = u_time; 
    float d30 = sdBettaTailShape3edges(rotate30 * center30, c30, r30, itime30);
   
   vec2 center31 = centcoord - vec2(0.230,0.450);
    float theta31 = 0.512;
    mat2 rotate31 = mat2(cos(theta31), -sin(theta31), sin(theta31), cos(theta31));
    float he31 = 0.304;
    float ra31 = 0.092;
    float rb31 = 0.224;
    float bu31 = 0.344;
    float d31 = sdEgg(rotate31 * center31, he31, ra31, rb31, bu31);

       float d3 = smin_cubic(d31, d30, 0.148);

    

    // --- 4. COMBINE ---
    float smooth_k = 0.096;
    
    float d_headbody = smoothUnion(d0, d1, smooth_k);
    float d_headbodyfin = smoothUnion(d_headbody, d3, 0.04);
    float d_tail = d2;
    float d_temp = smoothUnion(d_headbodyfin, d_tail, smooth_k);

    float mask = step(d_temp,0.0);
    
    // Output color
    vec3 colorwoutAlpha = vec3(centcoord * 0.5 + 0.5, 1.0);
    colour_out = mask * vec4(colorwoutAlpha, 1.0);
}