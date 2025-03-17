#include <raylib.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define WIN_W 600
#define WIN_H 800

const float wall_thickness = 20;
const float pad_standoff = 40;
const float ball_radius = 5;

const float speed_limit = 800;

#define COLNUM 14
#define ROWNUM 39

const int row_max = ROWNUM - 1;
const int col_max = COLNUM - 1;

const float block_width = (WIN_W - 2*wall_thickness) / COLNUM;
const float block_height = (WIN_H - wall_thickness) / ROWNUM;

struct box {
    float x;
    float y;
    float width;
    float height;
};

struct box left_wall = {wall_thickness / 2, WIN_H / 2 + wall_thickness, wall_thickness, WIN_H - wall_thickness};
struct box right_wall = {WIN_W - wall_thickness / 2, WIN_H / 2 + wall_thickness, wall_thickness, WIN_H - wall_thickness};
struct box ceiling = {WIN_W / 2, wall_thickness / 2, WIN_W - 2*wall_thickness, wall_thickness};

struct pad {
    float x;
    float y;
    float length;
    float thickness;
    float xmotion;
    int ghostly;
};

struct pad pad = {WIN_W/2, WIN_H - pad_standoff, 100, 20, 0, 0};

struct ball {
    float x;
    float y;
    float vx;
    float vy;
    int enable;
    int crushed;
};

//struct ball ball = {50,WIN_H-pad_standoff,0,0,1};
struct ball ball = {50,400,100,-400,1};

struct block {
    int enable : 1;
    int color : 4;
};

struct block block_matrix[ROWNUM][COLNUM];

void audio_callback(void *buffer, unsigned int frames);
void setup(void);
void mainloop_body(void);

void audio_callback(void *buffer, unsigned int frames)
{
    short *out = (short *)buffer;

    for (unsigned int i = 0; i < frames; i++)
    {
        out[i] = 0;
    }
}

float randf(void) {
    return (float) rand() / RAND_MAX;
}

Vector2 reflect(Vector2 N, Vector2 v) {
    float L = sqrt(N.x*N.x + N.y*N.y);
    if(L == 0) {
        printf("can't reflect using zero vector\n");
        abort();
    }
    float NNx = N.x / L;
    float NNy = N.y / L;
    float factor = 2 * (NNx * v.x + NNy * v.y);
    float dx = NNx * factor;
    float dy = NNy * factor;
    Vector2 ans = {v.x - dx, v.y - dy};
    return ans;
}

Vector2 clamp_length(Vector2 v, float limit) {
    float L = sqrt(v.x*v.x + v.y*v.y);
    if(L > limit){
        Vector2 vlimited = {v.x / L * limit, v.y / L * limit};
        printf("clamp_length triggered (speed limit), new value = %f %f\n", vlimited.x, vlimited.y);
        return vlimited;
    }
    else {
        return v;
    }
}

void draw_arena(void) {
    Color c = DARKGRAY;
    float thick = wall_thickness;
    DrawRectangle(0, 0, WIN_W, thick, c);
    DrawRectangle(0, 0, thick, WIN_H, c);
    DrawRectangle(WIN_W - thick, 0, thick, WIN_H, c);
}

void draw_pad(void) {
    float x = pad.x;
    float y = pad.y;
    float l = pad.length;
    float thick = pad.thickness;

    DrawRectangle(x - l/2, y - thick/2, l, thick, RAYWHITE);
}

void draw_ball(void) {
    DrawCircle(ball.x, ball.y, ball_radius, RAYWHITE);
}

void draw_box(struct box *bl, Color c) {
    DrawRectangle(bl->x - bl->width/2, bl->y - bl->height/2, bl->width, bl->height, c);
}

void draw_box_lines(struct box *bl, Color c) {
    DrawRectangleLines(bl->x - bl->width/2, bl->y - bl->height/2, bl->width, bl->height, c);
}

void draw_grid(){
    const int vspace = 20;
    const int hspace = 40;
    int rows = WIN_H / vspace;
    int cols = WIN_W / hspace;
    for(int i = 0; i < rows; i++) {
        DrawLine(0, i*vspace + wall_thickness, WIN_W, i*vspace + wall_thickness, DARKGRAY);
    }
    for(int j = 0; j < cols; j++) {
        DrawLine(j*hspace + wall_thickness, 0, j*hspace + wall_thickness, WIN_H, DARKGRAY);
    }
}



void player_control(void) {
    Vector2 motion = GetMouseDelta();
    pad.x += motion.x;
    pad.xmotion = motion.x;

    float left = wall_thickness;
    float right = WIN_W - wall_thickness;
    float half = pad.length / 2;

    if(pad.x < left + half) pad.x = left + half;
    if(pad.x > right - half) pad.x = right - half;
}

Vector2 normal_field_of_box(struct box b, float world_x, float world_y) {
    float x = world_x - b.x;
    float y = world_y - b.y;
    float left  = -b.width/2;
    float right = b.width/2;
    float top   = -b.height/2;
    float bottom = b.height/2;

    Vector2 N = {x, y};

    if(left < x && x < right && top < y && y < bottom) {
        N.x = 0;
        N.y = 0;
    }
    else if(left < x && x < right) {
        N.x = 0;
        N.y /= fabs(N.y);
    }
    else if(top < y && y < bottom) {
        N.y = 0;
        N.x /= fabs(N.x);
    }
    else {
        float aspect = b.width / b.height;
        N.y *= aspect;
        float L = sqrt(N.x*N.x + N.y*N.y);
        N.x /= L;
        N.y /= L;
    }

    return N;
}

int test_ball_box(struct ball * ball, struct box * box) {
    Vector2 center = {ball->x, ball->y};
    float radius = ball_radius;
    struct Rectangle rect = {box->x - box->width/2, box->y - box->height/2, box->width, box->height};
    return CheckCollisionCircleRec(center, radius, rect);
}

int test_ball_edge(struct ball * ball, Vector2 p0, Vector2 p1) {
    Vector2 center = {ball->x, ball->y};
    return CheckCollisionCircleLine(center, ball_radius, p0, p1);
}

/* assuming ball overlaps box, return a normal vector */
Vector2 test_ball_edges(struct ball * ball, struct box * box) {
    float left = box->x - box->width/2;
    float right = box->x + box->width/2;
    float top = box->y - box->height/2;
    float bottom = box->y + box->height/2;
    Vector2 a = {left, top};
    Vector2 b = {right, top};
    Vector2 c = {right, bottom};
    Vector2 d = {left, bottom};
    Vector2 normal = {0,0};
    if(test_ball_edge(ball, a, b)){ normal.y -= 1; }
    if(test_ball_edge(ball, b, c)){ normal.x += 1; }
    if(test_ball_edge(ball, c, d)){ normal.y += 1; }
    if(test_ball_edge(ball, d, a)){ normal.x -= 1; }
    return normal;
}

int ball_crush_test(int x, int y) {
    struct box pad_box = {pad.x, pad.y, pad.length, pad.thickness};
    struct ball ball = {x, y, 0, 0, 1};

    int touch_pad = 0;
    int touch_left = 0;
    int touch_right = 0;

    if(pad.ghostly == 0 && test_ball_box(&ball, &pad_box)) touch_pad = 1;
    if(test_ball_box(&ball, &left_wall)) touch_left = 1;
    if(test_ball_box(&ball, &right_wall)) touch_right = 1;

    return touch_pad && (touch_left || touch_right);
}

Vector2 quantize_velocity(Vector2 v) {
    float L = sqrt(v.x * v.x + v.y * v.y);
    float angle = atan2(v.y, v.x);

    if(0 <= angle && angle < PI/6) angle = PI/6;
    if(angle > 5*PI/6) angle = 5*PI/6;
    if(-PI/6 < angle && angle <= 0) angle = -PI/6;
    if(angle < -5*PI/6) angle = -5*PI/6;

    float sign = angle < 0 ? -1 : 1;
    float alpha = fabs(angle) / PI * 180;

    alpha = round((alpha - 30) / 10) * 10 + 30;
    alpha *= sign;
    angle = alpha/180 * PI;

    L = round(L / 100) * 100;
    if(L < 100) L = 100;
    if(L > 1000) L = 1000;

    Vector2 ans = {L * cos(angle), L * sin(angle)};
    return ans;
}


struct box generate_box_at(int row, int col) {
    float left = wall_thickness + col * block_width;
    float top  = wall_thickness + row * block_height;
    float x = left + block_width/2;
    float y = top + block_height/2;
    struct box b = {x, y, block_width, block_height};
    return b;
}

void neighborhood(float x, float y, Vector2 *norm_out, int *num_hits, Vector2 *blocks_hit, int *num_blocks, int *hit_pad) {

    struct ball ball = {x, y, 0, 0};

    *num_blocks = 0;
    *num_hits = 0;

    int row = floor((y - wall_thickness) / block_height);
    int col = floor((x - wall_thickness) / block_width);
    if(row < 0) return;
    if(col < 0) return;
    if(row > row_max) return;
    if(col > col_max) return;

    int bi = 0;
    Vector2 normal = {0,0};

    for(int i = -1; i <= 1; i++) {
        for(int j = -1; j <= 1; j++) {
            int a = row + i;
            int b = col + j;
            if(a < 0 || a > row_max) continue;
            if(b < 0 || b > col_max) continue;

            struct box box = generate_box_at(a, b);
            //draw_box_lines(&box, YELLOW);

            if(block_matrix[a][b].enable){
                if(test_ball_box(&ball, &box)) {
                    Vector2 n = test_ball_edges(&ball, &box);
                    normal.x += n.x;
                    normal.y += n.y;
                    (*num_blocks) += 1;
                    (*num_hits) += 1;
                    if(*num_hits > 4){
                        abort();
                    }

                    blocks_hit[bi].x = a;
                    blocks_hit[bi].y = b;
                    bi++;
                }
            }
        }
    }

    if(test_ball_box(&ball, &ceiling)){
        Vector2 n = test_ball_edges(&ball, &ceiling);
        normal.x += n.x;
        normal.y += n.y;
        (*num_hits) += 1;
    }
        
    if(test_ball_box(&ball, &left_wall)){
        Vector2 n = test_ball_edges(&ball, &left_wall);
        normal.x += n.x;
        normal.y += n.y;
        (*num_hits) += 1;
    }

    if(test_ball_box(&ball, &right_wall)){
        Vector2 n = test_ball_edges(&ball, &right_wall);
        normal.x += n.x;
        normal.y += n.y;
        (*num_hits) += 1;
    }

    struct box pad_box = {pad.x, pad.y, pad.length, pad.thickness};

    if(pad.ghostly == 0 && test_ball_box(&ball, &pad_box)){
        Vector2 n = test_ball_edges(&ball, &pad_box);
        normal.x += n.x;
        normal.y += n.y;
        (*num_hits) += 1;
        *hit_pad = 1;
    }

    norm_out->x = normal.x;
    norm_out->y = normal.y;

}

void ball_physics(void) {
    if(ball.enable == 0) return;
    if(ball.crushed) {
        ball.y += 1; 
        return;
    }

    if(ball.y > WIN_H) ball.y = 400;
    if(ball.x < 0) ball.x = 300;
    if(ball.x > WIN_W) ball.x = 300;

    float time_limit = 1.0 / 60;
    float x = ball.x;
    float y = ball.y;
    float vx = ball.vx;
    float vy = ball.vy;
    float speed = sqrt(vx*vx + vy*vy);

    while(time_limit > 0) {
        float time_chunk = 1.0 / speed;

        if(time_chunk > time_limit) time_chunk = time_limit;

        x += time_chunk * vx;
        y += time_chunk * vy;

        if(ball_crush_test(x, y)) {
            ball.crushed = 1;
            return;
        }

        int num_hits = 0;
        int num_blocks = 0;
        int hit_pad = 0;
        Vector2 normal = {0,0};
        Vector2 blocks_hit[4];
        neighborhood(x, y, &normal, &num_hits, blocks_hit, &num_blocks, &hit_pad);

        if(num_hits > 0) {
            //printf("COLLISION n = %d, normal = %f %f\n", num_blocks, normal.x, normal.y);
            if(normal.x == 0 && normal.y == 0) {
                printf("hit something but there's no normal\n");
                abort();
            }
            x += normal.x;
            y += normal.y;

            if(hit_pad) {
                Vector2 virtual_v = {vx, vy};
                Vector2 rv = reflect(normal, virtual_v);
                float bump = pad.xmotion * 10;
                int limit = 100;
                if(bump > limit) bump = limit;
                if(bump < -limit) bump = -limit;
//                printf("hit pad, xmotion = %f, bump = %f\n", pad.xmotion, bump);
                vx = rv.x + bump;
                /* correction and bounding pane */
                vy = rv.y;
            }
            else {
                float mag0 = sqrt(vx*vx + vy*vy);
                Vector2 rv = reflect(normal, (Vector2){vx,vy});
                vx = rv.x;
                vy = rv.y;
                float mag1 = sqrt(vx*vx + vy*vy);
                if(fabs(mag0 - mag1) > 0.1) {
                    printf("magnitude changed = %f => %f\n", mag0, mag1);
                    abort();
                }
            }


            if(num_hits > 4) abort();
            for(int i = 0; i < num_blocks; i++) {
                block_matrix[(int)blocks_hit[i].x][(int)blocks_hit[i].y].enable = 0;
            }
        }

        time_limit -= time_chunk;
    }

    ball.x = x;
    ball.y = y;
    ball.vx = vx;
    ball.vy = vy;

    Vector2 qv = quantize_velocity((Vector2){ball.vx, ball.vy});
    ball.vx = qv.x;
    ball.vy = qv.y;

}


// move the pad by delta_x very carefully
void pad_physics(float delta_x) {
    if(isnan(delta_x)){
        printf("delta_x is NaN\n");
        abort();
    }

    delta_x = trunc(delta_x);

    if(delta_x == 0) return;

    float increment = delta_x < 0 ? -1 : 1;

    //printf("pad physics delta_x = %f, increment = %f\n", delta_x, increment);



    int watchdog = 0;

    while(delta_x != 0) {
        //printf("delta_x = %f\n", delta_x);
        pad.x += increment;
        delta_x -= increment;

        struct box pad_box = {pad.x, pad.y, pad.length, pad.thickness};

        if(pad.ghostly == 0 && test_ball_box(&ball, &pad_box)){
            Vector2 normal = test_ball_edges(&ball, &pad_box);

            ball.x += normal.x;
            ball.y += normal.y;

            Vector2 virtual_v = {ball.vx - pad.xmotion, ball.vy};

            Vector2 rv = reflect(normal, virtual_v);
            ball.vx = rv.x + pad.xmotion;
            ball.vy = rv.y;

        }
        
        watchdog++;
        if(watchdog > 10000) {
            printf("pad_physics watchdog triggered. delta_x = %f, increment = %f\n", delta_x, increment);
            abort();
        }

    }

    if(pad.x - pad.length/2 < wall_thickness) pad.x = wall_thickness + pad.length / 2;
    if(pad.x + pad.length/2 > WIN_W - wall_thickness) pad.x = WIN_W - wall_thickness - pad.length / 2;
}



void setup_blocks() {

    for(int i = 0; i < ROWNUM; i++) {
        for(int j = 0; j < COLNUM; j++) {
            block_matrix[i][j].enable = 0;
            block_matrix[i][j].color = 0;
        }
    }

    for(int i = 0; i < 15; i++) {
        for(int j = 0; j < COLNUM; j++) {
            int n = 2*i + j;
            block_matrix[i][j].color = n % 7 + 1;
            block_matrix[i][j].enable = 1;
        }
    }

}

void draw_blocks(){
    for(int i = 0; i < ROWNUM; i++) {
        for(int j = 0; j < COLNUM; j++) {
            if(block_matrix[i][j].enable == 0) continue;
            struct box box = generate_box_at(i,j);
            Color c = GRAY;
            switch(block_matrix[i][j].color){
                case 0: c = DARKGRAY; break;
                case 1: c = RED; break;
                case 2: c = GREEN; break;
                case 3: c = BLUE; break;
                case 4: c = SKYBLUE; break;
                case 5: c = PURPLE; break;
                case 6: c = DARKGREEN; break;
                case 7: c = ORANGE; break;
            }
            draw_box(&box, c);
        }
    }
    
}

void setup() {
    setup_blocks();
}

void mainloop_body(void) {
    Vector2 mouse = GetMousePosition();

    // CTRL
    //player_control();

    // FIZIKS
    ball_physics();

    if(IsKeyPressed(KEY_G)){
        ball.vx *= 2;
        ball.vy *= 2;
    }

    Vector2 delta = GetMouseDelta();
    pad.xmotion = delta.x;
    pad_physics(delta.x);

    // GRAPHICS
    BeginDrawing();

    ClearBackground(BLACK);

    //draw_grid();

    draw_arena();
    draw_pad();
    draw_ball();

    float speed = sqrt(pow(ball.vx, 2) + pow(ball.vy,2));
    //DrawText(TextFormat("vx = %f", ball.vx), 30, WIN_H - 100, 20, WHITE);
    //DrawText(TextFormat("vy = %f", ball.vy), 30, WIN_H -  80, 20, WHITE);
    //DrawText(TextFormat("v  = %f", speed),   30, WIN_H -  60, 20, WHITE);

    if(ball.crushed) {
        DrawText("BALL CRUSHED", 150, 600, 40, WHITE);
    }

/*
    struct box test = generate_box_at(2, 10);
    Vector2 N = normal_field_of_box(test, ball.x, ball.y);
    DrawLine(ball.x, ball.y, ball.x + 30*N.x, ball.y + 30*N.y, YELLOW);
    draw_box(test, SKYBLUE);


    draw_box(generate_box_at(2, 3), GREEN);
    draw_box(generate_box_at(3, 3), YELLOW);
    draw_box(generate_box_at(4, 5), RED);
    draw_box(generate_box_at(5, 6), BLUE);
*/

    draw_blocks();

    EndDrawing();
}

int main(int argc, char * argv[]) {

    InitWindow(WIN_W, WIN_H, "BREAKOUT");

    DisableCursor();
    SetTargetFPS(60);

    InitAudioDevice();
    SetAudioStreamBufferSizeDefault(4096);
    AudioStream stream = LoadAudioStream(44100, 16, 1);
    SetAudioStreamCallback(stream, audio_callback);
    PlayAudioStream(stream);

    setup();

    while(WindowShouldClose() == 0) mainloop_body();

    UnloadAudioStream(stream);
    CloseAudioDevice();

    CloseWindow();

    return 0;

}
