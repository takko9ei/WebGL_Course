//this file to backup implemented glsl code
//need too place these code into a real proj

// credit: jikihakase (will show up on website)
// butterfly: Monarch Butterfly
// credit: jikihakase (will show up on website)
// butterfly: Monarch Butterfly
uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;

float sdEgg( in vec2 p, in float he, in float ra, in float rb, in float bu );
float sdCircle( vec2 p, float r );

float sdOrientedVesica( vec2 p, vec2 a, vec2 b, float w );
float smoothUnion(float d1, float d2, float k);
float sdBettaTailShape( vec2 p, vec2 c, float r, float iTime );

// 辅助函数：简单的噪声组合，用于生成不规则的褶皱波形
// angle: 当前的角度
// time: 时间变量，让褶皱自身轻微蠕动


void main() {
   vec2 centcoord = gl_FragCoord.xy/u_resolution.xy;
   
   float theta3 = -4.360;
   mat2 rotate3 = mat2(cos(theta3),-sin(theta3),sin(theta3),cos(theta3));
   vec2 center3 = centcoord - vec2(0.570,0.570);
   vec2 c3 = vec2(0.920,0.280);
   float r3 = 0.500;
   float itime3 = 0.200;
   float d3 = sdBettaTailShape(rotate3*center3,c3,r3,itime3);
   
   
   vec2 center2 = centcoord-vec2(0.650,0.460);
   vec2 a2 = vec2(0.320,0.380);
   vec2 b2 = vec2(0.070,0.060);
   float w2 = 0.088;
   float d2 =sdOrientedVesica(center2,a2,b2,w2);
   
   //**rule here: the no. of param of shape #n, follows its no
   //egg shape 0
   vec2 center0 = centcoord-vec2(0.820,0.670);//center definition, same below
   float theta0 = -4.256;
   mat2 rotate0 = mat2(cos(theta0),-sin(theta0),sin(theta0),cos(theta0));
   float he0 = 0.244; 
    // ra: Bottom radius (radius at the base)
    float ra0 = 0.102;
    // rb: Top radius (radius at the tip)
    float rb0 = 0.072; 
    // bu: "Bulge" factor? (Usually 1.0 for standard geometry calculation)
    //     Normally this variable is not in standard sdEgg, but based on formula, 
    //     it scales 'r'. Keep it at 1.0 or adjust slightly for shape variations.
    float bu0 = 0.264;
   float d0 = sdEgg(rotate0*center0, he0, ra0, rb0, bu0);
   
   vec2 center1 = centcoord-vec2(0.520,0.740);
   float radius1 = 0.236;
   float d1 = sdCircle(center1,radius1);
   
   //phase blendshape
   
   
   float d = d3;
   d = smoothUnion(d,d2,0.048);
   d = smoothUnion(d,d0,0.048);
   
   //here, d1 no use
   //implemented the body part
   //the body part will be one texture
   
   //and the fin part will be one
   //the tail part will be one
   //total 3 parts
   
   //optional: fore fin
   
   
    float mask = step(d, 0.0);
   colour_out = vec4( centcoord,0.0, 1.0 );
   colour_out *= mask;
}
float getRuffleNoise(float angle, float time) {
    // 叠加三层正弦波
    // 层1：主要的大波浪 (低频高幅)
    float noise = sin(angle * 8.0 + time * 1.5) * 0.5;
    // 层2：中等褶皱 (中频中幅)
    noise += sin(angle * 17.0 - time * 2.0) * 0.25;
    // 层3：细碎的边缘纹理 (高频低幅)
    noise += sin(angle * 37.0 + time * 3.5) * 0.125;
    
    // 将结果归一化到大致 [-1, 1] 区间，然后乘以一个强度系数
    return noise * 0.1; // 0.1 控制整个褶皱的深度
}

// 修改后的扇形 SDF
float sdBettaTailShape( vec2 p, vec2 c, float r, float iTime )
{
    // --- 新增部分：计算荷叶边偏移量 ---
    // 1. 在 abs(p.x) 之前获取原始角度。atan(y, x) 返回 [-PI, PI]
    float angle = atan(p.y, p.x);

    // 2. 计算基于角度的半径偏移量。
    // 我们乘以 r 是为了让褶皱的大小与尾巴整体大小成比例。
    float ruffleOffset = getRuffleNoise(angle, iTime) * r;

    // --- 原始 sdPie 逻辑开始 ---
    p.x = abs(p.x);
    
    // 3. 【关键修改】将固定的半径 r 替换为动态的半径 (r + ruffleOffset)
    // 这样只有外侧的弧形边缘会变得不规则，而两侧的直边保持不变。
    float l = length(p) - (r + ruffleOffset);
    
    float m = length(p-c*clamp(dot(p,c),0.0,r)); // c=sin/cos of aperture
    return max(l,m*sign(c.y*p.x-c.x*p.y));
}

float sdEgg( in vec2 p, in float he, in float ra, in float rb, in float bu )
{
    // all this can be precomputed for any given shape
    float r = 0.5*(he + ra+rb)/bu;
    float da = r - ra;
    float db = r - rb;
    float y = (db*db - da*da - he*he)/(2.0*he);
    float x = sqrt(da*da - y*y);
     
    // only this needs to be run per pixel
    p.x = abs(p.x);
    float k = p.y*x - p.x*y;
    if( k>0.0 && k<he*(p.x+x) )
        return length(p+vec2(x,y))-r;
    return min( length(p)-ra,
                length(vec2(p.x,p.y-he))-rb );
}

float sdCircle( vec2 p, float r )
{
    return length(p) - r;
}

//helpers
//smooth boool
float smoothUnion(float d1, float d2, float k) {//k represents the factor of smooth
    // Calculate interpolation factor 'h' based on the difference between d1 and d2.
    // clamp() ensures 'h' stays within the valid range [0.0, 1.0].
    float h = clamp(0.5 + 0.5 * (d2 - d1) / k, 0.0, 1.0);
    
    // Linearly interpolate between d2 and d1 using 'h'.
    // Then subtract a correction term 'k * h * (1.0 - h)' to smooth out the sharp corner.
    return mix(d2, d1, h) - k * h * (1.0 - h);
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

//float intersect
//float except