#include "raylib.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace
{
    std::size_t idx(int x, int y, int n)
    {
        return static_cast<std::size_t>(x + y * n);
    }

    void set_bnd(int b, std::vector<float> &x, int n)
    {
        for (int i = 1; i < n - 1; ++i)
        {
            x[idx(0, i, n)] = b == 1 ? -x[idx(1, i, n)] : x[idx(1, i, n)];
            x[idx(n - 1, i, n)] = b == 1 ? -x[idx(n - 2, i, n)] : x[idx(n - 2, i, n)];
            x[idx(i, 0, n)] = b == 2 ? -x[idx(i, 1, n)] : x[idx(i, 1, n)];
            x[idx(i, n - 1, n)] = b == 2 ? -x[idx(i, n - 2, n)] : x[idx(i, n - 2, n)];
        }

        x[idx(0, 0, n)] = 0.5f * (x[idx(1, 0, n)] + x[idx(0, 1, n)]);
        x[idx(0, n - 1, n)] = 0.5f * (x[idx(1, n - 1, n)] + x[idx(0, n - 2, n)]);
        x[idx(n - 1, 0, n)] = 0.5f * (x[idx(n - 2, 0, n)] + x[idx(n - 1, 1, n)]);
        x[idx(n - 1, n - 1, n)] = 0.5f * (x[idx(n - 2, n - 1, n)] + x[idx(n - 1, n - 2, n)]);
    }

    void lin_solve(int b, std::vector<float> &x, const std::vector<float> &x0, float a, float c, int iter, int n)
    {
        const float c_recip = 1.0f / c;

        for (int k = 0; k < iter; ++k)
        {
            for (int y = 1; y < n - 1; ++y)
            {
                for (int x_pos = 1; x_pos < n - 1; ++x_pos)
                {
                    x[idx(x_pos, y, n)] =
                        (x0[idx(x_pos, y, n)] + a * (x[idx(x_pos + 1, y, n)] + x[idx(x_pos - 1, y, n)] + x[idx(x_pos, y + 1, n)] + x[idx(x_pos, y - 1, n)])) * c_recip;
                }
            }

            set_bnd(b, x, n);
        }
    }

    void diffuse(int b, std::vector<float> &x, const std::vector<float> &x0, float diff, float dt, int iter, int n)
    {
        const float a = dt * diff * static_cast<float>((n - 2) * (n - 2));
        lin_solve(b, x, x0, a, 1.0f + 4.0f * a, iter, n);
    }

    void project(std::vector<float> &veloc_x, std::vector<float> &veloc_y, std::vector<float> &p, std::vector<float> &div, int iter, int n)
    {
        for (int y = 1; y < n - 1; ++y)
        {
            for (int x_pos = 1; x_pos < n - 1; ++x_pos)
            {
                div[idx(x_pos, y, n)] = -0.5f * (veloc_x[idx(x_pos + 1, y, n)] - veloc_x[idx(x_pos - 1, y, n)] + veloc_y[idx(x_pos, y + 1, n)] - veloc_y[idx(x_pos, y - 1, n)]) / static_cast<float>(n);
                p[idx(x_pos, y, n)] = 0.0f;
            }
        }

        set_bnd(0, div, n);
        set_bnd(0, p, n);
        lin_solve(0, p, div, 1.0f, 4.0f, iter, n);

        for (int y = 1; y < n - 1; ++y)
        {
            for (int x_pos = 1; x_pos < n - 1; ++x_pos)
            {
                veloc_x[idx(x_pos, y, n)] -= 0.5f * (p[idx(x_pos + 1, y, n)] - p[idx(x_pos - 1, y, n)]) * static_cast<float>(n);
                veloc_y[idx(x_pos, y, n)] -= 0.5f * (p[idx(x_pos, y + 1, n)] - p[idx(x_pos, y - 1, n)]) * static_cast<float>(n);
            }
        }

        set_bnd(1, veloc_x, n);
        set_bnd(2, veloc_y, n);
    }

    void advect(int b, std::vector<float> &d, const std::vector<float> &d0, const std::vector<float> &veloc_x, const std::vector<float> &veloc_y, float dt, int n)
    {
        const float dtx = dt * static_cast<float>(n - 2);
        const float dty = dt * static_cast<float>(n - 2);
        const float max_coord = static_cast<float>(n) - 1.5f;

        for (int y = 1; y < n - 1; ++y)
        {
            for (int x_pos = 1; x_pos < n - 1; ++x_pos)
            {
                float x = static_cast<float>(x_pos) - dtx * veloc_x[idx(x_pos, y, n)];
                float y_pos = static_cast<float>(y) - dty * veloc_y[idx(x_pos, y, n)];

                x = std::clamp(x, 0.5f, max_coord);
                y_pos = std::clamp(y_pos, 0.5f, max_coord);

                const int x0 = static_cast<int>(std::floor(x));
                const int x1 = x0 + 1;
                const int y0 = static_cast<int>(std::floor(y_pos));
                const int y1 = y0 + 1;

                const float s1 = x - static_cast<float>(x0);
                const float s0 = 1.0f - s1;
                const float t1 = y_pos - static_cast<float>(y0);
                const float t0 = 1.0f - t1;

                d[idx(x_pos, y, n)] =
                    s0 * (t0 * d0[idx(x0, y0, n)] + t1 * d0[idx(x0, y1, n)]) +
                    s1 * (t0 * d0[idx(x1, y0, n)] + t1 * d0[idx(x1, y1, n)]);
            }
        }

        set_bnd(b, d, n);
    }
}

struct Cell
{
    int size{};
    float dt{};
    float diff{};
    float visc{};

    std::vector<float> s;
    std::vector<float> density;
    std::vector<float> vx;
    std::vector<float> vy;
    std::vector<float> vx0;
    std::vector<float> vy0;

    Cell(int gridSize, float timeStep, float diffusion, float viscosity)
        : size(gridSize),
          dt(timeStep),
          diff(diffusion),
          visc(viscosity),
          s(static_cast<std::size_t>(gridSize * gridSize), 0.0f),
          density(static_cast<std::size_t>(gridSize * gridSize), 0.0f),
          vx(static_cast<std::size_t>(gridSize * gridSize), 0.0f),
          vy(static_cast<std::size_t>(gridSize * gridSize), 0.0f),
          vx0(static_cast<std::size_t>(gridSize * gridSize), 0.0f),
          vy0(static_cast<std::size_t>(gridSize * gridSize), 0.0f)
    {
    }

    void addDensity(int x, int y, float amount)
    {
        density[idx(x, y, size)] += amount;
    }

    void addVelocity(int x, int y, float amount_x, float amount_y)
    {
        const std::size_t index = idx(x, y, size);
        vx[index] += amount_x;
        vy[index] += amount_y;
    }

    void step()
    {
        diffuse(1, vx0, vx, visc, dt, 4, size);
        diffuse(2, vy0, vy, visc, dt, 4, size);

        project(vx0, vy0, vx, vy, 4, size);

        advect(1, vx, vx0, vx0, vy0, dt, size);
        advect(2, vy, vy0, vx0, vy0, dt, size);

        project(vx, vy, vx0, vy0, 4, size);

        diffuse(0, s, density, diff, dt, 4, size);
        advect(0, density, s, vx, vy, dt, size);
    }
};

int main()
{
    const int N = 200;          // N x N grid
    const int aspect_ratio = 4; // each cell is aspect_ratio x aspect_ratio pixels

    Cell cell(N, 0.1f, 0.0001f, 0.0001f);

    InitWindow(N * aspect_ratio, N * aspect_ratio, "Stam99");
    SetTargetFPS(30);

    // Graphics
    Texture2D density_texture = LoadTextureFromImage(GenImageColor(N, N, BLACK));
    SetTextureFilter(density_texture, TEXTURE_FILTER_POINT); // prevent blurring when scaling up

    std::vector<Color> pixels(static_cast<std::size_t>(N * N));

    const int source_band_height = 3;
    const int source_band_width = N / 10;
    const int source_band_start_x = (N - source_band_width) / 2;
    const int source_band_end_x = source_band_start_x + source_band_width;
    const float source_density = 0.5f;
    const float source_velocity_y = -2.0f;

    while (!WindowShouldClose())
    {
        for (int y = N - source_band_height; y < N - 1; ++y)
        {
            for (int x = source_band_start_x; x < source_band_end_x; ++x)
            {
                cell.addDensity(x, y, source_density);
                cell.addVelocity(x, y, 0.0f, source_velocity_y);
            }
        }

        cell.step();

        for (int i = 0; i < N * N; ++i)
        {
            float d = std::clamp(cell.density[i], 0.0f, 1.0f);
            unsigned char a = static_cast<unsigned char>(255.0f * d);
            pixels[static_cast<std::size_t>(i)] = Color{0, 0, 0, a};
        }

        UpdateTexture(density_texture, pixels.data());

        BeginDrawing();
        ClearBackground(RAYWHITE);

        DrawTexturePro(
            density_texture,
            Rectangle{0.0f, 0.0f, static_cast<float>(N), static_cast<float>(N)},
            Rectangle{0.0f, 0.0f, static_cast<float>(N * aspect_ratio), static_cast<float>(N * aspect_ratio)},
            Vector2{0.0f, 0.0f},
            0.0f,
            WHITE);

        DrawFPS(10, 10);
        EndDrawing();
    }

    CloseWindow();
    return 0;
}