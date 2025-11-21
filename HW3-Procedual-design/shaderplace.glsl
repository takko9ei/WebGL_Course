//this file to backup implemented glsl code
//need too place these code into a real proj

// credit: jikihakase (will show up on website)
// butterfly: Monarch Butterfly
uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;

float sdEgg( in vec2 p, in float he, in float ra, in float rb, in float bu );
float sdCircle( vec2 p, float r );
float sdOrientedVesica( vec2 p, vec2 a, vec2 b, float w );
float smoothUnion(float d1, float d2, float k);

void main() {
   vec2 centcoord = gl_FragCoord.xy/u_resolution.xy;
   
   vec2 center2 = centcoord-vec2(0.580,0.540);
   vec2 a2 = vec2(0.320,0.380);
   vec2 b2 = vec2(0.070,0.060);
   float w2 = 0.088;
   float d2 =sdOrientedVesica(center2,a2,b2,w2);
   
   //**rule here: the no. of param of shape #n, follows its no
   //egg shape 0
   vec2 center0 = centcoord-vec2(0.730,0.730);//center definition, same below
   float theta0 = -4.256;
   mat2 rotate0 = mat2(cos(theta0),-sin(theta0),sin(theta0),cos(theta0));
   float he0 = 0.236; 
    // ra: Bottom radius (radius at the base)
    float ra0 = 0.102;
    // rb: Top radius (radius at the tip)
    float rb0 = 0.056; 
    // bu: "Bulge" factor? (Usually 1.0 for standard geometry calculation)
    //     Normally this variable is not in standard sdEgg, but based on formula, 
    //     it scales 'r'. Keep it at 1.0 or adjust slightly for shape variations.
    float bu0 = 0.296;
   float d0 = sdEgg(rotate0*center0, he0, ra0, rb0, bu0);
   
   vec2 center1 = centcoord-vec2(0.520,0.740);
   float radius1 = 0.236;
   float d1 = sdCircle(center1,radius1);
   
   //phase blendshape
   
   
   float d = smoothUnion(d0,d2,0.048);
   
   
    float mask = step(d, 0.0);
   colour_out = vec4( centcoord,0.0, 1.0 );
   colour_out *= mask;
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