# oop

This directory serves as a quick reference for object-oriented programming in R.

While everything is an object in R, many of the objects commonly used in R is behave quite differently than objects in general purpose programming languages, such as C++, Java and Python.

R have 3 object oriented systems: S3, S4, and Reference Classes (aka S5).

S3 objects are very informal, any object can be given a new class by:
d = list(name="Daisy", age="14", breed="shepard")
class(d) = "dog"

While the ad hoc nature of the S3 class system can be very useful, it can also make standarization of similar ojects very difficult. Instead, it can be helpful to use constructor functions to create objects:

dog = function(n,a,b) {
value = list(name=n, age=a, breed=b)
attr(value, "class") = "dog"
value
}
