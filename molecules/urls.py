from django.urls import path
from . import views

urlpatterns = [
    path('', views.search_molecule, name='search_molecule'),
]
