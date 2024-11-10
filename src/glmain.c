#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <stdio.h>
#include <math.h>

// Window dimensions
const unsigned int WIDTH = 800;
const unsigned int HEIGHT = 600;

// Player position and angle
float playerX = 0.0f, playerY = 1.0f, playerZ = 5.0f;
float playerAngle = 0.0f;

// Movement and rotation speed
float moveSpeed = 0.1f;
float rotationSpeed = 1.5f;

// Ground size for demo
float groundSize = 20.0f;

// Function prototypes
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window);
void renderScene();
void drawCube(float size);
void setPerspectiveProjection(float fov, float aspect, float near, float far);

// Initialize OpenGL settings
void initGL() {
    glEnable(GL_DEPTH_TEST); // Enable depth testing
    glClearColor(0.5f, 0.7f, 0.9f, 1.0f); // Set sky-blue background
}

// Function to update camera position
void updateCamera() {
    glRotatef(-playerAngle, 0.0f, 1.0f, 0.0f);
    glTranslatef(-playerX, -playerY, -playerZ);
}

// Render the ground and cube
void renderScene() {
    // Draw ground
    glColor3f(0.3f, 0.7f, 0.3f); // Green ground
    glBegin(GL_QUADS);
    glVertex3f(-groundSize, 0.0f, -groundSize);
    glVertex3f(groundSize, 0.0f, -groundSize);
    glVertex3f(groundSize, 0.0f, groundSize);
    glVertex3f(-groundSize, 0.0f, groundSize);
    glEnd();

    // Draw a simple cube as a landmark
    glPushMatrix();
    glTranslatef(0.0f, 0.5f, -5.0f);
    glColor3f(1.0f, 0.0f, 0.0f); // Red cube
    drawCube(1.0f);
    glPopMatrix();
}

// Draws a simple cube
void drawCube(float size) {
    float halfSize = size / 2.0f;
    glBegin(GL_QUADS);

    // Front face
    glVertex3f(-halfSize, -halfSize, halfSize);
    glVertex3f(halfSize, -halfSize, halfSize);
    glVertex3f(halfSize, halfSize, halfSize);
    glVertex3f(-halfSize, halfSize, halfSize);

    // Back face
    glVertex3f(-halfSize, -halfSize, -halfSize);
    glVertex3f(halfSize, -halfSize, -halfSize);
    glVertex3f(halfSize, halfSize, -halfSize);
    glVertex3f(-halfSize, halfSize, -halfSize);

    // Top face
    glVertex3f(-halfSize, halfSize, -halfSize);
    glVertex3f(halfSize, halfSize, -halfSize);
    glVertex3f(halfSize, halfSize, halfSize);
    glVertex3f(-halfSize, halfSize, halfSize);

    // Bottom face
    glVertex3f(-halfSize, -halfSize, -halfSize);
    glVertex3f(halfSize, -halfSize, -halfSize);
    glVertex3f(halfSize, -halfSize, halfSize);
    glVertex3f(-halfSize, -halfSize, halfSize);

    // Left face
    glVertex3f(-halfSize, -halfSize, -halfSize);
    glVertex3f(-halfSize, halfSize, -halfSize);
    glVertex3f(-halfSize, halfSize, halfSize);
    glVertex3f(-halfSize, -halfSize, halfSize);

    // Right face
    glVertex3f(halfSize, -halfSize, -halfSize);
    glVertex3f(halfSize, halfSize, -halfSize);
    glVertex3f(halfSize, halfSize, halfSize);
    glVertex3f(halfSize, -halfSize, halfSize);

    glEnd();
}

// Set a perspective projection matrix
void setPerspectiveProjection(float fov, float aspect, float near, float far) {
    float top = tan(fov / 2.0f * M_PI / 180.0f) * near;
    float right = top * aspect;

    glFrustum(-right, right, -top, top, near, far);
}

int main() {
    // Initialize GLFW
    if (!glfwInit()) {
        fprintf(stderr, "Failed to initialize GLFW\n");
        return -1;
    }

    // Set OpenGL version to 3.3 Core
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // Create window
    GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "3D Game with GLFW", NULL, NULL);
    if (!window) {
        fprintf(stderr, "Failed to create GLFW window\n");
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // Load OpenGL using GLAD
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        fprintf(stderr, "Failed to initialize GLAD\n");
        return -1;
    }

    // Initialize OpenGL settings
    initGL();

    // Set the perspective projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    setPerspectiveProjection(60.0f, (float)WIDTH / (float)HEIGHT, 0.1f, 100.0f);
    glMatrixMode(GL_MODELVIEW);

    // Main loop
    while (!glfwWindowShouldClose(window)) {
        // Input
        processInput(window);

        // Render
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glLoadIdentity();
        updateCamera();
        renderScene();

        // Swap buffers and poll events
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // Clean up and exit
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}

// Process keyboard input
void processInput(GLFWwindow* window) {
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        playerX += sin(playerAngle * M_PI / 180.0f) * moveSpeed;
        playerZ -= cos(playerAngle * M_PI / 180.0f) * moveSpeed;
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        playerX -= sin(playerAngle * M_PI / 180.0f) * moveSpeed;
        playerZ += cos(playerAngle * M_PI / 180.0f) * moveSpeed;
    }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        playerAngle -= rotationSpeed;
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        playerAngle += rotationSpeed;
    }
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
        glfwSetWindowShouldClose(window, 1);
    }
}

// Adjust viewport on window resize
void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    setPerspectiveProjection(60.0f, (float)width / (float)height, 0.1f, 100.0f);
    glMatrixMode(GL_MODELVIEW);
}
