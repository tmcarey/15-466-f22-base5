#include "Scene.hpp"

struct Ghost {
	Scene::Transform *transform;
	bool isFound = false;
	glm::vec3 startPosition;
};
