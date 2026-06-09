#include "raylib.h"

int main()
{
    const int screenWidth = 800;
    const int screenHeight = 450;

    InitWindow(screenWidth, screenHeight, "Lava Lamp");

    Vector2 ballPosition = {GetScreenWidth() / 2.0f, GetScreenHeight() / 2.0f};
    Vector2 ballSpeed = {0.0f, 4.0f};
    int ballRadius = 20;
    float gravity = 0.2f;

    bool useGravity = true;
    bool pause = 0;
    int framesCounter = 0;

    SetTargetFPS(60);

    while (!WindowShouldClose())
    {
        ballPosition.x += ballSpeed.x;
        ballPosition.y += ballSpeed.y;

        if (useGravity)
            ballSpeed.y += gravity;

        if ((ballPosition.x >= (GetScreenWidth() - ballRadius)) || (ballPosition.x <= ballRadius))
            ballSpeed.x *= -1.0f;
        if ((ballPosition.y >= (GetScreenHeight() - ballRadius)) || (ballPosition.y <= ballRadius))
            ballSpeed.y *= -0.95f;

        BeginDrawing();

        ClearBackground(RAYWHITE);

        DrawCircleV(ballPosition, (float)ballRadius, MAROON);

        DrawFPS(10, 10);

        EndDrawing();
    }

    CloseWindow();

    return 0;
}