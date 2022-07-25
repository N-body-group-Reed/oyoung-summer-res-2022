#include <glad/glad.h> /* must add glad.c */
#include <GLFW/glfw3.h>
#include <glm/glm/glm.hpp>
#include <glm/glm/gtc/matrix_transform.hpp>
#include <glm/glm/gtc/type_ptr.hpp>

#include <iostream>
#include <string>
#include <cmath>
#include <vector>

#include "Shader.h"
#include "Texture.h"
#include "Camera.h"
#include "Geometry.h"

#include <NSim.h> /* Must still add cpp files for library */

const double pi = 3.14159265358979323846;

// function declaration
void glfwInitSetup(void);
GLFWwindow* glfwWindowSetup(void);
void gladInit(void);

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);



// NSim graphics functions
void printOCTree(class OCTree* T, class Shader& S, unsigned int bBox, unsigned int p, glm::mat4& projection, glm::vec4& color);
void print_octree(class OCTree* T, class Shader& S, unsigned int bBox, unsigned int p,glm::vec3 center,int depth,int index);
//void printTree(class OCTree* T, class Shader& S, unsigned int bBox, unsigned int p, glm::mat4& projection, glm::vec4& color);
//void print_tree(class OCTree* T, class Shader& S, unsigned int bBox, unsigned int p,glm::vec3 center,int depth,int index);




// constants
#define SCR_WIDTH 1000
#define SCR_HEIGHT 800

// global
float alpha_v = 0.0f;
// camera
Camera C(glm::vec3(0.0f, 0.0f, 15.0f));
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;
// timing
float deltaTime;
float lastFrame;

int main() {
    glfwInitSetup();

    GLFWwindow* window = glfwWindowSetup();

    gladInit();

    glEnable(GL_DEPTH_TEST);

    int nrAttributes;
    glGetIntegerv(GL_MAX_VERTEX_ATTRIBS, &nrAttributes);
    std::cout << "Maximum vertex attributes supported: " << nrAttributes << std::endl;

    Shader S1("VertexShader.vs", "FragmentShader.fs");
    Shader S2("VertexSimple.vs", "FragmentSimple.fs");

    Sphere S(.5, 40, 40);

    std::vector<float> vs = S.vertices;
    std::vector<unsigned int> indices = S.indices;

    unsigned int S_VBO, S_VAO, S_EBO;
    glGenVertexArrays(1, &S_VAO);
    glGenBuffers(1, &S_VBO);
    glGenBuffers(1, &S_EBO);

    glBindVertexArray(S_VAO);

    glBindBuffer(GL_ARRAY_BUFFER, S_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vs[0]) * vs.size(), vs.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, S_EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices[0]) * indices.size(), indices.data(), GL_STATIC_DRAW);

    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    // draw in wireframe polygons.
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    RectangularPrism RP(1, 1, 1, Decomp_type::WIR);

    vs = RP.vertices;
    indices = RP.indices;

    unsigned int RP_VBO, RP_VAO, RP_EBO;
    glGenVertexArrays(1, &RP_VAO);
    glGenBuffers(1, &RP_VBO);
    glGenBuffers(1, &RP_EBO);

    glBindVertexArray(RP_VAO);

    glBindBuffer(GL_ARRAY_BUFFER, RP_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vs[0]) * vs.size(), vs.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, RP_EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices[0]) * indices.size(), indices.data(), GL_STATIC_DRAW);

    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    // draw in wireframe polygons.
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    Rectangle R(1, 1, Decomp_type::WIR);

    vs = R.vertices;
    indices = R.indices;

    unsigned int R_VBO, R_VAO, R_EBO;
    glGenVertexArrays(1, &R_VAO);
    glGenBuffers(1, &R_VBO);
    glGenBuffers(1, &R_EBO);

    glBindVertexArray(R_VAO);

    glBindBuffer(GL_ARRAY_BUFFER, R_VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vs[0]) * vs.size(), vs.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, R_EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices[0]) * indices.size(), indices.data(), GL_STATIC_DRAW);
    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    // draw in wireframe polygons.
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glm::vec4 orange(1.0f, 0.5f, 0.2f, 1.0f);
    glm::vec4 purple(1.0f, 0.0f, 0.5f, 1.0f);
    glm::vec4 green(0.2f, 0.7f, 0.4f, 1.0f);





    
    // NSim
    std::vector<class Particle> ps_tree, ps_direct;
    NSim_Init(ps_tree, 1000);
    NSim_Init(ps_direct, 1000);
    std::list<class Particle*> ptrs;
    std::list<class Particle*>::iterator left, right, part;
    double dt = 1/100.0;
    double before;
    double after;




    
    lastFrame = static_cast<float>(glfwGetTime());
    while (!glfwWindowShouldClose(window)) {
        // per-frame time logic
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
	//std::cout<<"       TIME"<<deltaTime<<std::endl;
        lastFrame = currentFrame;

        //std::cout << 1/deltaTime << std::endl;

        processInput(window);

        // render background
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glm::mat4 projection = glm::perspective(glm::radians(C.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        
        S2.set();
        S2.setUniform_Mat4("projection", projection);
        S2.setUniform_Mat4("view", C.GetViewMatrix());







	
        // NSim
	/*

for (std::vector<class Particle>::iterator p = ps_tree.begin(); p != ps_tree.end(); p++) ptrs.push_back(&*p);
	before = static_cast<float>(glfwGetTime());
        left = ptrs.begin();
        right = ptrs.end();
        Tree* T = new Tree(ptrs, left, right,0,0,"octree");
        part = T->Partition(ptrs,left,right,0,0);
	left = ptrs.begin();
        right = ptrs.end();
	
        buildTree(T, ptrs, left, right, part);
        T->computeMassMoments(left, right);
        orderParticles(ps_tree, ptrs);

        //printTree(T, S2, RP_VAO, R_VAO, projection, orange);

        NSim_Step(ps_tree, T, dt);
        delete T;	
	after = static_cast<float>(glfwGetTime());
	ptrs.erase(ptrs.begin(), ptrs.end());
std::cout<<"         OCTIME "<<after-before<<std::endl;
	//*/


	
        
	for (std::vector<class Particle>::iterator p = ps_tree.begin(); p != ps_tree.end(); p++) ptrs.push_back(&*p);
	before = static_cast<float>(glfwGetTime());
        left = ptrs.begin();
        right = ptrs.end();
        Tree* T = new Tree(ptrs, left, right,0,0,"octree");
        part = T->getPartitionIterator(left,right);	
        buildTree(T, ptrs, left, right, part);
        T->computeMassMoments(left, right);
        orderParticles(ps_tree, ptrs);

        //printTree(T, S2, RP_VAO, R_VAO, projection, orange);

        NSim_Step(ps_tree, T, dt);
        delete T;	
	after = static_cast<float>(glfwGetTime());
	ptrs.erase(ptrs.begin(), ptrs.end());
	std::cout<<"         TREETIME "<<after-before<<std::endl;
	//*/

	
	for (std::vector<class Particle>::iterator p = ps_direct.begin(); p != ps_direct.end(); p++) ptrs.push_back(&*p);
	before = static_cast<float>(glfwGetTime());
	OCTree* D = new OCTree(ptrs, left, right,0,0,"direct sum");
        NSim_Step(ps_direct, D,dt);
	delete D;
	after = static_cast<float>(glfwGetTime());
	std::cout<<"         TIME "<<after-before<<std::endl;
	//*/
	ptrs.erase(ptrs.begin(), ptrs.end());



	
	
	
        



	
        glBindVertexArray(S_VAO);
        S2.setUniform_Vec4("color", purple);
        for (int i = 0; i < ps_tree.size(); i++) {
            glm::vec3 offset(ps_tree[i].coords[0], ps_tree[i].coords[1], ps_tree[i].coords[2]);
            glm::mat4 model = glm::mat4(1.0f);
            model = glm::translate(model, offset);
            S2.setUniform_Mat4("model", model);
            glDrawElements(GL_TRIANGLES, S.indices.size(), GL_UNSIGNED_INT, 0);
        }
        S2.setUniform_Vec4("color", green);
        for (int i = 0; i < ps_direct.size(); i++) {
            glm::vec3 offset(ps_direct[i].coords[0], ps_direct[i].coords[1], ps_direct[i].coords[2]);
            glm::mat4 model = glm::mat4(1.0f);
            model = glm::translate(model, offset);
            S2.setUniform_Mat4("model", model);
            glDrawElements(GL_TRIANGLES, S.indices.size(), GL_UNSIGNED_INT, 0);
	    }

        glfwSwapBuffers(window);
        glfwPollEvents();
	
    }

    //glDeleteVertexArrays(1, &VAO);
    //glDeleteBuffers(1, &VBO);

    glDeleteVertexArrays(1, &RP_VAO);
    glDeleteBuffers(1, &RP_VBO);
    glDeleteBuffers(1, &RP_EBO);

    glDeleteVertexArrays(1, &R_VAO);
    glDeleteBuffers(1, &R_VBO);
    glDeleteBuffers(1, &R_EBO);

    glDeleteVertexArrays(1, &S_VAO);
    glDeleteBuffers(1, &S_VBO);
    glDeleteBuffers(1, &S_EBO);

    glfwTerminate();
    return 0;
}

// function implementations

void glfwInitSetup(void) {
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); // 4
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3); // 6
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
}

GLFWwindow* glfwWindowSetup(void) {
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Test", NULL, NULL);
    if (window == NULL) {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return NULL;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    return window;
}

void gladInit(void) {
    // glad: load all OpenGL function pointers
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << std::endl;
        glfwTerminate();
    }
}

void processInput(GLFWwindow* window) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
    else if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS) {
        alpha_v += .1f;
        if (alpha_v >= 1.0f) alpha_v = 1.0f;
    }
    else if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS) {
        alpha_v -= .1f;
        if (alpha_v <= 0.0f) alpha_v = 0.0f;
    }

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        C.ProcessKeyboard(Camera_Movement::FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        C.ProcessKeyboard(Camera_Movement::BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        C.ProcessKeyboard(Camera_Movement::LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        C.ProcessKeyboard(Camera_Movement::RIGHT, deltaTime);
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}

// glfw: whenever the mouse moves, this callback is called
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn) {
    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);

    if (firstMouse) {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top
    lastX = xpos;
    lastY = ypos;

    C.ProcessMouseMovement(xoffset, yoffset);
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    C.ProcessMouseScroll(static_cast<float>(yoffset));
}










// NSim

void printOCTree(class OCTree* T, class Shader& S, unsigned int bBox, unsigned int p, glm::mat4& projection, glm::vec4& color) {
    // enable shader
    S.set();
    S.setUniform_Mat4("projection", projection);
    S.setUniform_Mat4("view", C.GetViewMatrix());
    S.setUniform_Vec4("color", color);
    print_octree(T, S, bBox, p, glm::vec3(0.0,0.0,12.0),0,0);
}

void print_octree(class OCTree* T, class Shader& S, unsigned int bBox, unsigned int p,glm::vec3 center,int depth,int index) {
    // check if leaf contains only one particle
  if(index>2){
    index=0;
    depth+=1;
  }
  std::cout<<index<<","<<depth<<std::endl;
  glm::mat4 model = glm::mat4(1.0f);
  glm::vec3 s =glm::vec3 (pow(0.5,depth),pow(0.5,depth),pow(0.5,depth));
  glm::vec3 t;
  if (index==0){t=glm::vec3(pow(0.5,depth+2),0,0);
  }else if (index==1) {t=glm::vec3(0,pow(0.5,depth+2),0);
  }else {t=glm::vec3(0,0,pow(0.5,depth+2));}
  glBindVertexArray(bBox);
  model = glm::translate(model, glm::vec3(center));
  model = glm::scale(model, glm::vec3(s));
  S.setUniform_Mat4("model", model);
  if(index==0)glDrawElements(GL_LINES, 24, GL_UNSIGNED_INT, 0);
  // recurse
  if (T->l != nullptr) print_octree(T->l, S, bBox, p,center+t,depth,index+1);
  if (T->r != nullptr) print_octree(T->r, S, bBox, p,center,depth,index+1);
}
