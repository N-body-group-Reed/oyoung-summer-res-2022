/*

void printTree(class Tree* T, class Shader& S, unsigned int bBox, unsigned int p, glm::mat4& projection, glm::vec4& color) {
    // enable shader
    S.set();
    S.setUniform_Mat4("projection", projection);
    S.setUniform_Mat4("view", C.GetViewMatrix());
    S.setUniform_Vec4("color", color);
    print_tree(T, S, bBox, p);
}

void print_tree(class Tree* T, class Shader& S, unsigned int bBox, unsigned int p) {
    // check if leaf contains only one particle
    double tol = 1e-5;
    glm::mat4 model = glm::mat4(1.0f);
    if (T->x_max - T->x_min > tol && T->y_max - T->y_min > tol && T->z_max - T->z_min > tol) {
        // draw bounding box
        glBindVertexArray(bBox);
        model = glm::translate(model, glm::vec3((T->x_max + T->x_min) / 2, (T->y_max + T->y_min) / 2, (T->z_max + T->z_min) / 2));
        model = glm::scale(model, glm::vec3(T->x_max - T->x_min, T->y_max - T->y_min, T->z_max - T->z_min));
        S.setUniform_Mat4("model", model);
        glDrawElements(GL_LINES, 24, GL_UNSIGNED_INT, 0);
    }
    if (T->l != nullptr || T->r != nullptr) {
        // draw partition
        glBindVertexArray(p);
        model = glm::mat4(1.0f);
        // determine axis
        if (T->axis == Partition_axis::X) {
            model = glm::rotate(model, 3 * (float)pi / 2, glm::vec3(0, 1, 0));
            model = glm::translate(model, glm::vec3((T->z_max + T->z_min) / 2, (T->y_max + T->y_min) / 2, -(T->x_max + T->x_min) / 2));
            model = glm::scale(model, glm::vec3(T->z_max - T->z_min, T->y_max - T->y_min, 1));
        }
        else if (T->axis == Partition_axis::Y) {
            model = glm::rotate(model, (float)pi / 2, glm::vec3(1, 0, 0));
            model = glm::translate(model, glm::vec3((T->x_max + T->x_min) / 2, (T->z_max + T->z_min) / 2, -(T->y_max + T->y_min) / 2));
            model = glm::scale(model, glm::vec3(T->x_max - T->x_min, T->z_max - T->z_min, 1));
        }
        else {
            model = glm::translate(model, glm::vec3((T->x_max + T->x_min) / 2, (T->y_max + T->y_min) / 2, (T->z_max + T->z_min) / 2));
            model = glm::scale(model, glm::vec3(T->x_max - T->x_min, T->y_max - T->y_min, 1));
        }
        S.setUniform_Mat4("model", model);
        glDrawElements(GL_LINE_LOOP, 4, GL_UNSIGNED_INT, 0); // glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
    }
    // recurse
    if (T->l != nullptr) print_tree(T->l, S, bBox, p);
    if (T->r != nullptr) print_tree(T->r, S, bBox, p);
}

*/